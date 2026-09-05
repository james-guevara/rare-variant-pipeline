data "aws_caller_identity" "current" {}

resource "aws_ecr_repository" "targeted" {
  name                 = var.ecr_repository_name
  image_tag_mutability = "IMMUTABLE"

  image_scanning_configuration {
    scan_on_push = true
  }
}

resource "aws_ecr_lifecycle_policy" "targeted" {
  repository = aws_ecr_repository.targeted.name
  policy = jsonencode({
    rules = [{
      rulePriority = 1
      description  = "Retain the newest 30 commit images"
      selection = {
        tagStatus     = "tagged"
        tagPrefixList = ["sha-"]
        countType     = "imageCountMoreThan"
        countNumber   = 30
      }
      action = { type = "expire" }
    }]
  })
}

resource "aws_cloudwatch_log_group" "codebuild" {
  name              = "/aws/codebuild/${var.codebuild_project_name}"
  retention_in_days = 30
}

resource "aws_iam_role" "codebuild" {
  name = "${var.codebuild_project_name}-service"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "codebuild.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy" "codebuild" {
  role = aws_iam_role.codebuild.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "logs:CreateLogStream",
          "logs:PutLogEvents"
        ]
        Resource = "${aws_cloudwatch_log_group.codebuild.arn}:*"
      },
      {
        Effect   = "Allow"
        Action   = ["ecr:GetAuthorizationToken"]
        Resource = "*"
      },
      {
        Effect = "Allow"
        Action = [
          "ecr:BatchCheckLayerAvailability",
          "ecr:CompleteLayerUpload",
          "ecr:DescribeImages",
          "ecr:GetDownloadUrlForLayer",
          "ecr:InitiateLayerUpload",
          "ecr:PutImage",
          "ecr:UploadLayerPart"
        ]
        Resource = aws_ecr_repository.targeted.arn
      }
    ]
  })
}

resource "aws_codebuild_project" "targeted" {
  name         = var.codebuild_project_name
  service_role = aws_iam_role.codebuild.arn

  artifacts { type = "NO_ARTIFACTS" }
  source {
    type = "NO_SOURCE"
    buildspec = yamlencode({
      version = 0.2
      env = {
        "exported-variables" = ["IMAGE_URI", "IMAGE_DIGEST"]
        shell                = "bash"
      }
      phases = {
        build = {
          commands = [
            "git clone --no-checkout https://github.com/$GITHUB_REPOSITORY.git .",
            "git fetch --depth 1 origin $GITHUB_SHA",
            "git checkout --detach $GITHUB_SHA",
            "source scripts/ci/codebuild_targeted_container.sh"
          ]
        }
      }
    })
  }

  environment {
    compute_type                = var.codebuild_compute_type
    image                       = "aws/codebuild/standard:7.0"
    type                        = "LINUX_CONTAINER"
    image_pull_credentials_type = "CODEBUILD"
    privileged_mode             = true

    environment_variable {
      name  = "ECR_REPOSITORY"
      value = aws_ecr_repository.targeted.name
    }
  }

  logs_config {
    cloudwatch_logs {
      group_name  = aws_cloudwatch_log_group.codebuild.name
      stream_name = "container-build"
    }
  }

  build_timeout = 60
}

resource "aws_iam_openid_connect_provider" "github" {
  url            = "https://token.actions.githubusercontent.com"
  client_id_list = ["sts.amazonaws.com"]
}

resource "aws_iam_role" "github_actions" {
  name = "${var.codebuild_project_name}-github-actions"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect = "Allow"
      Principal = {
        Federated = aws_iam_openid_connect_provider.github.arn
      }
      Action = "sts:AssumeRoleWithWebIdentity"
      Condition = {
        StringEquals = {
          "token.actions.githubusercontent.com:aud" = "sts.amazonaws.com"
        }
        StringLike = {
          "token.actions.githubusercontent.com:sub" = "repo:${var.github_repository}:*"
        }
      }
    }]
  })
}

resource "aws_iam_role_policy" "github_actions" {
  role = aws_iam_role.github_actions.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect = "Allow"
      Action = [
        "codebuild:BatchGetProjects",
        "codebuild:BatchGetBuilds",
        "codebuild:StartBuild"
      ]
      Resource = aws_codebuild_project.targeted.arn
    }]
  })
}
