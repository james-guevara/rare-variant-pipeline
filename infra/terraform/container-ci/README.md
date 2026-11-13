# Targeted container CI

This small Terraform stack creates an immutable ECR repository, one ordinary
CodeBuild project, and a narrowly scoped GitHub OIDC role. It deliberately does
not use a persistent runner, GitHub webhook, PAT, or CodeBuild source connection.

## Bootstrap

Run:

```bash
cp backend.hcl.example backend.hcl
# Set the existing Terraform state bucket in backend.hcl.
terraform init
terraform validate
terraform plan
terraform apply
```

The stack creates GitHub's account-wide OIDC identity provider. If another
Terraform stack later needs to own that shared provider, import or move this
resource rather than creating a duplicate.

Copy the outputs into these GitHub repository variables:

- `AWS_CODEBUILD_ROLE_ARN` <- `github_actions_role_arn`
- `AWS_CODEBUILD_PROJECT` <- `codebuild_project_name`
- `AWS_REGION` <- `us-east-1`

The workflow is manual initially. It can be added as a required pull-request
check only after a successful build and comparison with the existing GHCR image.
