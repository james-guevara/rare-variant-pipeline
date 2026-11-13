variable "aws_region" {
  description = "AWS region used for builds and the container registry."
  type        = string
  default     = "us-east-1"
}

variable "github_repository" {
  description = "GitHub repository allowed to start builds, as owner/name."
  type        = string
  default     = "james-guevara/rare-variant-pipeline"
}

variable "ecr_repository_name" {
  description = "ECR repository for targeted annotation images."
  type        = string
  default     = "rare-variant-pipeline-targeted"
}

variable "codebuild_project_name" {
  description = "Name of the heavyweight container build project."
  type        = string
  default     = "rare-variant-pipeline-targeted-container"
}

variable "codebuild_compute_type" {
  description = "CodeBuild compute size. LARGE is a conservative starting point."
  type        = string
  default     = "BUILD_GENERAL1_LARGE"
}

variable "tags" {
  description = "Tags applied to CI resources."
  type        = map(string)
  default = {
    Project   = "rare-variant-pipeline"
    ManagedBy = "Terraform"
    Component = "container-ci"
  }
}
