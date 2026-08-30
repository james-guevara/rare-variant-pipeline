# G2MH indel consequence adjudication

Eight normalized insertion alleles produced different selected-transcript coding
consequences between Ensembl VEP 115 and FastVEP in the comprehensive G2MH run.
They were adjudicated independently from either consequence label by reconstructing
the Ensembl 115 transcript CDS from the chromosome GFF3 and reference FASTA,
inserting the normalized alternate sequence at its genomic breakpoint, and
translating the resulting CDS.

| Gene | GRCh38 allele | Transcript | Inserted bases | Breakpoint | Translation result | Adjudication |
|---|---|---|---:|---|---|---|
| MSH6 | chr2:47806618 T>TTGAGAAGATGAATCAGTCACTACGATTATTTCGGTAACTAACTAACTATAATGGAATTATAACTAACTGACCTTAAGTTTCAAAG | ENST00000234420 | 85 | CDS 3968/4083 | stop 1361 -> 1335 | frameshift |
| DOCK3 | chr3:51381547 A>AGGGGAGCAGTGAGGGGCAAC | ENST00000266037 | 20 | CDS 6081/6093 | protein through the existing stop is unchanged | stop-retained repeat insertion |
| DOCK2 | chr5:169674423 T>TCCAAAATTGACTATGGCAACAAGTAACCTCTCTTTCCTCTGCAAAGAGGTTTTCTTCCC | ENST00000520908 | 59 | CDS 448/5493 | stop 1831 -> 158 | frameshift |
| SORCS1 | chr10:106577495 T>TTCAC | ENST00000263054 | 4 | CDS 3431/3507 | stop 1169 -> 1145 | frameshift |
| MTMR2 | chr11:95849862 C>CTGGGATACGGCCTCTTGATCTGAAGGATGCCACTCTCTTTAATT | ENST00000346299 | 44 | outside CDS, at exon boundary | unchanged | splice-region/coding-boundary, not frameshift |
| DCC | chr18:53499510 G>GGTAAA | ENST00000442544 | 5 | outside CDS, at exon boundary | unchanged | splice-region/coding-boundary, not frameshift |
| PRDM12 | chr9:130668313 G>GGT | ENST00000253008 | 2 | outside CDS, at exon boundary | unchanged | splice-region/coding-boundary, not frameshift |
| PRDM12 | chr9:130668313 G>GGTGT | ENST00000253008 | 4 | outside CDS, at exon boundary | unchanged | splice-region/coding-boundary, not frameshift |

The result supports the FastVEP-selected consequences for all eight alleles. In
particular, insertion length modulo three is insufficient for repeat-associated
terminal variants such as DOCK3: the altered transcript must be translated through
the existing stop. No additional pipeline correction is required, and the corrected
comprehensive G2MH total remains `lof_t2 = 1212`.

This adjudication concerns these exact alleles and Ensembl release 115 transcripts.
The standalone LOFTEE classifier remains responsible for evaluating only genuine
LoF consequence terms after transcript selection.
