output "codebuild_project_name" {
  value = aws_codebuild_project.targeted.name
}

output "ecr_repository_url" {
  value = aws_ecr_repository.targeted.repository_url
}

output "github_actions_role_arn" {
  value = aws_iam_role.github_actions.arn
}
