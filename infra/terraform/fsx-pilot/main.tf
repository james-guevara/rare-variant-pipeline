resource "aws_fsx_lustre_file_system" "pilot" {
  storage_capacity         = var.storage_capacity_gib
  storage_type             = "SSD"
  subnet_ids               = [var.subnet_id]
  security_group_ids       = var.security_group_ids
  deployment_type          = "SCRATCH_2"
  file_system_type_version = "2.15"
  data_compression_type    = "LZ4"

  tags = {
    Name = var.name
  }
}

resource "aws_fsx_data_repository_association" "pilot_inputs" {
  file_system_id                   = aws_fsx_lustre_file_system.pilot.id
  file_system_path                 = var.file_system_path
  data_repository_path             = var.s3_data_repository_path
  batch_import_meta_data_on_create = true
  imported_file_chunk_size         = 1024
}

resource "aws_launch_template" "vcz_batch" {
  name        = "rare-variant-vcz-fsx-imported"
  description = "terraform-bootstrap"

  user_data = base64encode(<<-USER_DATA
    MIME-Version: 1.0
    Content-Type: multipart/mixed; boundary="==VCZFSX=="

    --==VCZFSX==
    Content-Type: text/x-shellscript; charset="us-ascii"

    #!/bin/bash
    set -euxo pipefail
    dnf install -y lustre-client
    mkdir -p /fsx
    echo '${aws_fsx_lustre_file_system.pilot.dns_name}@tcp:/${aws_fsx_lustre_file_system.pilot.mount_name} /fsx lustre defaults,noatime,flock,_netdev 0 0' >> /etc/fstab
    for attempt in $(seq 1 30); do
      if mountpoint -q /fsx || mount /fsx; then
        break
      fi
      sleep 10
    done
    mountpoint -q /fsx

    --==VCZFSX==--
  USER_DATA
  )

  tags = {
    Name = "rare-variant-vcz-fsx"
  }
}

resource "aws_batch_compute_environment" "vcz" {
  name         = "rare-variant-vcz-fsx"
  type         = "MANAGED"
  state        = "ENABLED"
  service_role = "arn:aws:iam::${data.aws_caller_identity.current.account_id}:role/aws-service-role/batch.amazonaws.com/AWSServiceRoleForBatch"

  compute_resources {
    type                = "EC2"
    allocation_strategy = "BEST_FIT_PROGRESSIVE"
    min_vcpus           = 0
    desired_vcpus       = 0
    max_vcpus           = var.batch_max_vcpus
    instance_type       = var.batch_instance_types
    subnets             = [var.subnet_id]
    security_group_ids  = var.security_group_ids
    instance_role       = var.batch_instance_profile_arn

    launch_template {
      launch_template_id = aws_launch_template.vcz_batch.id
      version            = "$Latest"
    }

    ec2_configuration {
      image_type = "ECS_AL2023"
    }

    tags = {
      Name = "rare-variant-vcz-fsx"
    }
  }

  lifecycle {
    create_before_destroy = true
  }
}

resource "aws_batch_job_queue" "vcz" {
  name     = "rare-variant-vcz-fsx"
  state    = "ENABLED"
  priority = 10

  compute_environment_order {
    order               = 1
    compute_environment = aws_batch_compute_environment.vcz.arn
  }
}

resource "aws_batch_job_definition" "vcz_smoke" {
  name                  = "rare-variant-vcz-fsx-smoke"
  type                  = "container"
  platform_capabilities = ["EC2"]

  container_properties = jsonencode({
    image      = "public.ecr.aws/docker/library/python:3.12-slim"
    command    = ["true"]
    jobRoleArn = "arn:aws:iam::${data.aws_caller_identity.current.account_id}:role/ecsInstanceRole"
    resourceRequirements = [
      { type = "VCPU", value = "1" },
      { type = "MEMORY", value = "1024" }
    ]
    volumes = [{
      name = "fsx"
      host = { sourcePath = "/fsx" }
    }]
    mountPoints = [{
      sourceVolume  = "fsx"
      containerPath = "/fsx"
      readOnly      = false
    }]
  })
}

resource "aws_batch_job_definition" "vcz_conversion" {
  name                  = "rare-variant-vcz-conversion"
  type                  = "container"
  platform_capabilities = ["EC2"]

  container_properties = jsonencode({
    image      = var.vcz_batch_image
    command    = ["/usr/local/bin/run_chromosome.sh", "Ref::chromosome"]
    jobRoleArn = "arn:aws:iam::${data.aws_caller_identity.current.account_id}:role/ecsInstanceRole"
    resourceRequirements = [
      { type = "VCPU", value = "16" },
      { type = "MEMORY", value = "57344" }
    ]
    volumes = [{
      name = "fsx"
      host = { sourcePath = "/fsx" }
    }]
    mountPoints = [{
      sourceVolume  = "fsx"
      containerPath = "/fsx"
      readOnly      = false
    }]
  })

  parameters = {
    chromosome = "21"
  }
}

resource "aws_batch_job_definition" "vcz_conversion_32" {
  name                  = "rare-variant-vcz-conversion-32"
  type                  = "container"
  platform_capabilities = ["EC2"]

  container_properties = jsonencode({
    image      = var.vcz_batch_image
    command    = ["/usr/local/bin/run_chromosome.sh", "Ref::chromosome"]
    jobRoleArn = "arn:aws:iam::${data.aws_caller_identity.current.account_id}:role/ecsInstanceRole"
    environment = [
      { name = "VCZ_WORKERS", value = "32" },
      { name = "VCZ_MEMORY", value = "96G" },
    ]
    resourceRequirements = [
      { type = "VCPU", value = "32" },
      { type = "MEMORY", value = "114688" }
    ]
    volumes = [{
      name = "fsx"
      host = { sourcePath = "/fsx" }
    }]
    mountPoints = [{
      sourceVolume  = "fsx"
      containerPath = "/fsx"
      readOnly      = false
    }]
  })

  parameters = {
    chromosome = "1"
  }
}

data "aws_caller_identity" "current" {}
