# Disposable FSx pilot

This Terraform root creates the temporary filesystem and its matching compute layer:

1. An FSx for Lustre `SCRATCH_2` filesystem.
2. A data repository association (DRA) that exposes an existing S3 prefix inside FSx.
3. An EC2 launch template that mounts that exact filesystem at `/fsx`.
4. A 16-vCPU-capable AWS Batch compute environment and queue.
5. A small FSx mount smoke-test job definition.

The existing EC2 instance, IAM roles, subnet, security group, S3 objects, and state bucket are inputs.
They are not imported into this state and cannot be deleted by destroying this root.

## Why this shape

S3 is durable storage. FSx is a temporary, high-throughput working view of the pilot
prefix. Existing S3 object metadata is imported when the DRA is created; file contents
are fetched lazily on first access and then remain warm on Lustre for repeated tests.
No automatic export policy is configured, so benchmark outputs should be copied to S3
explicitly before the filesystem is destroyed.

## Run from the persistent EC2 host

Create private `backend.hcl` and `terraform.tfvars` files from the examples, then:

```bash
terraform init -backend-config=backend.hcl
terraform fmt -check -recursive
terraform validate
terraform plan -out=create.tfplan
terraform show create.tfplan
```

Apply only the reviewed saved plan:

```bash
terraform apply create.tfplan
terraform output -raw mount_command
```

Install the Lustre client, create `/fsx`, and run the emitted mount command. Confirm that
the linked path contains the expected S3 objects before hydrating the large pVCF.

## Cleanup boundary

Before destroy, confirm that no process uses `/fsx`, copy durable outputs to S3, and
compare object counts and sizes. Then create and review a saved destroy plan. Never treat
the FSx filesystem as the only copy of an input or result.
