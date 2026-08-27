#!/usr/bin/env perl
use strict;
use warnings;
use Storable qw(fd_retrieve);
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;

print join("\t", qw(
  transcript_id transcript_version mane_select mane_plus_clinical canonical
  appris tsl protein_coding ccds length ensembl refseq
)), "\n";

my %seen;
for my $path (@ARGV) {
  open my $input, '-|', 'gzip', '-dc', $path
    or die "cannot open $path: $!\n";
  my $cache = fd_retrieve($input);
  for my $seqname (keys %$cache) {
    for my $transcript (@{$cache->{$seqname} || []}) {
      my $id = $transcript->stable_id;
      next unless $id && !$seen{$id}++;
      my %attribute;
      for my $attr (@{$transcript->get_all_Attributes}) {
        push @{$attribute{$attr->code}}, $attr->value;
      }
      my $mane = exists $attribute{MANE_Select} ? 1 : 0;
      my $mane_plus = exists $attribute{MANE_Plus_Clinical} ? 1 : 0;
      my $tsl = 100;
      if (exists $attribute{TSL} && $attribute{TSL}[0] =~ /tsl(\d+)/) {
        $tsl = $1;
      }
      my $appris = 100;
      if (exists $attribute{appris} && $attribute{appris}[0] =~ /([A-Za-z]).+(\d+)/) {
        my ($type, $grade) = ($1, $2);
        $grade += 10 if lc(substr($type, 0, 1)) eq 'a';
        $appris = $grade;
      }
      my $length;
      if ($transcript->translation) {
        $length = eval {
          length(
            $transcript->{_variation_effect_feature_cache}->{translateable_seq}
              || $transcript->translateable_seq
          )
        };
      }
      $length = $transcript->length unless defined $length;
      my $source = lc($transcript->{_source_cache} || '');
      print join("\t",
        $id,
        $transcript->version || '',
        $mane,
        $mane_plus,
        $transcript->is_canonical ? 1 : 0,
        $appris,
        $tsl,
        $transcript->biotype eq 'protein_coding' ? 1 : 0,
        ($transcript->{_ccds} && $transcript->{_ccds} ne '-') ? 1 : 0,
        0 - $length,
        $source eq 'ensembl' ? 0 : 1,
        $source eq 'refseq' ? 0 : 1,
      ), "\n";
    }
  }
}
