# Software environment and execution provenance

## Paragraph genotyping

The production workflow used Paragraph 2.3
from Bioconda.

Captured Paragraph environment:
- Python 3.12.3
- pysam 0.22.1
- One Paragraph thread per sample

The 21 newly processed samples used individually
calibrated maximum-read settings (-M) between
486 and 803.

The 803m sample reused the earlier Paragraph
pilot batches containing 72, 6 and 18 graphs.

## Reference

Genome: Capra hircus ARS1.2
Assembly: GCF_001704415.2

RefSeq annotation: Release 102

GFF3 SHA-256:
b702d44f0ebdd638100bc7b93a08163a4e35a61e6651cb398e49d542d5a0b826

## Validation

22 animals x 96 loci = 2,112 genotype attempts.

FT-PASS calls: 1,645
Conservative loci: 47
Variable conservative loci: 39

Independent annotation reconstruction:
96/96 classifications matched.

The original annotation source script was not
recovered. The independent reannotation is a
validation, not a claim of historical script identity.

## Reproducibility limitations

Input BAMs, FASTQs, reference genome files and
individual-level source genotype data are not
included in this software repository.

The three historical production scripts are
preserved in src/ with verified SHA-256 hashes.
