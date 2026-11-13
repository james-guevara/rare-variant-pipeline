variable "aws_region" {
  description = "AWS region containing the pilot EC2 instance and S3 bucket."
  type        = string
  default     = "us-east-1"
}

variable "name" {
  description = "Name tag for the disposable filesystem."
  type        = string
  default     = "rare-variant-targeting-pilot"
}

variable "subnet_id" {
  description = "Existing subnet in the same Availability Zone as the EC2 client."
  type        = string
}

variable "security_group_ids" {
  description = "Existing security groups that allow Lustre traffic between clients and FSx."
  type        = list(string)
}

variable "storage_capacity_gib" {
  description = "FSx SCRATCH_2 capacity. AWS permits 1200, 2400, then increments of 2400 GiB."
  type        = number
  default     = 1200

  validation {
    condition = (
      var.storage_capacity_gib == 1200 ||
      var.storage_capacity_gib == 2400 ||
      (var.storage_capacity_gib > 2400 && var.storage_capacity_gib % 2400 == 0)
    )
    error_message = "Use 1200, 2400, or a larger multiple of 2400 GiB."
  }
}

variable "s3_data_repository_path" {
  description = "Existing S3 bucket or prefix exposed through the Lustre namespace."
  type        = string

  validation {
    condition     = startswith(var.s3_data_repository_path, "s3://")
    error_message = "s3_data_repository_path must start with s3://."
  }
}

variable "file_system_path" {
  description = "Absolute directory inside FSx associated with the S3 prefix."
  type        = string
  default     = "/rare-variant-pilot"

  validation {
    condition     = startswith(var.file_system_path, "/")
    error_message = "file_system_path must be absolute."
  }
}

variable "batch_instance_profile_arn" {
  description = "IAM instance profile used by EC2 hosts in the FSx-enabled Batch environment."
  type        = string
}

variable "batch_max_vcpus" {
  description = "Maximum aggregate vCPUs for the FSx-enabled Batch compute environment."
  type        = number
  default     = 64
}

variable "vcz_batch_image" {
  description = "Immutable ECR image used by VCZ conversion jobs."
  type        = string
}

variable "targeted_annotation_image" {
  description = "Immutable ECR image digest used by the targeted FastVEP/LOFTEE workflow."
  type        = string
}

variable "batch_instance_types" {
  description = "EC2 instance types allowed for 16-, 32-, and 64-vCPU VCZ conversion jobs."
  type        = list(string)
  default = [
    "m7i.4xlarge", "m6i.4xlarge", "m7a.4xlarge", "m6a.4xlarge",
    "m7i.8xlarge", "m6i.8xlarge", "m7a.8xlarge", "m6a.8xlarge",
    "m7i.16xlarge", "m6i.16xlarge", "m7a.16xlarge", "m6a.16xlarge",
  ]
}

variable "tags" {
  description = "Tags applied to resources managed by this pilot."
  type        = map(string)
  default = {
    Project   = "rare-variant-pipeline"
    ManagedBy = "Terraform"
    Lifecycle = "temporary"
  }
}
