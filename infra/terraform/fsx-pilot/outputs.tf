output "file_system_id" {
  description = "ID of the disposable FSx for Lustre filesystem."
  value       = aws_fsx_lustre_file_system.pilot.id
}

output "dns_name" {
  description = "DNS name used by Lustre clients."
  value       = aws_fsx_lustre_file_system.pilot.dns_name
}

output "mount_name" {
  description = "Lustre mount name assigned by FSx."
  value       = aws_fsx_lustre_file_system.pilot.mount_name
}

output "mount_command" {
  description = "Command to mount this filesystem after the Lustre client is installed."
  value       = "sudo mount -t lustre -o relatime,flock ${aws_fsx_lustre_file_system.pilot.dns_name}@tcp:/${aws_fsx_lustre_file_system.pilot.mount_name} /fsx"
}

output "linked_input_path" {
  description = "Path under the /fsx mount backed by the configured S3 prefix."
  value       = "/fsx${var.file_system_path}"
}

output "data_repository_association_id" {
  description = "ID of the S3 data repository association."
  value       = aws_fsx_data_repository_association.pilot_inputs.id
}

output "vcz_batch_queue" {
  description = "FSx-enabled AWS Batch queue for VCZ conversion jobs."
  value       = aws_batch_job_queue.vcz.name
}

output "vcz_batch_smoke_job_definition" {
  description = "Job definition used to verify that Batch hosts can see the FSx mount."
  value       = aws_batch_job_definition.vcz_smoke.arn
}

output "vcz_batch_conversion_job_definition" {
  description = "Pinned job definition for chromosome VCZ conversion."
  value       = aws_batch_job_definition.vcz_conversion.arn
}

output "vcz_batch_conversion_32_job_definition" {
  description = "Thirty-two-worker VCZ conversion job definition ARN."
  value       = aws_batch_job_definition.vcz_conversion_32.arn
}

output "targeted_chr22_job_definition" {
  description = "Fork-pinned FastVEP/LOFTEE chr22 job definition ARN."
  value       = aws_batch_job_definition.targeted_chr22.arn
}
