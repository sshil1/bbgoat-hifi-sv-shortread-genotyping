# Black Bengal goat HiFi-informed structural-variant genotyping (22 Illumina genomes)

This public repository accompanies an **original-research manuscript in preparation**, extending a separate two-PacBio-HiFi-animal SV discovery Data Note. The study targeted 96 previously discovered alleles (48 DEL and 48 INS, stratified by size) using Paragraph 2.3 in 22 Black Bengal goat short-read BAMs. **This repository's existing GitHub V0.2.0 release and Zenodo record are historical snapshots; the main branch contains later documentation corrections.**

## Independently audited findings

| Measure | Result |
|---|---:|
| Genotype attempts | 2,112 (22 animals × 96 SVs) |
| FT-PASS calls | 1,645 |
| FT-PASS calls with DP ≥ 10 | 1,632 |
| Loci with ≥ 20 passing calls | 52 |
| Conservative analysis loci | 47 |
| Variable genotype-state loci among 47 | 39 |
| Sensitivity analysis excluding SV0078 | 46 loci; 39 variable |

The selected panel is **ascertained from two HiFi genomes**. FT-PASS is a callability screen, not genotype truth or clinical validation.

## Methods and provenance

- 21 Illumina samples have archived executed Paragraph commands, manifests and calibration JSONs; the 803m sample reused three earlier 72+6+18-graph pilot batches. A separate four-graph test was a smoke test, not four additional study loci.
- Executed Paragraph commands used one thread and calibrated read caps of 486–803. The 21 runs completed 96 graphs apiece; 11 read-limit warnings appeared in six sample logs.
- The original source-selection script and production controller exist in the authors' verified server-side software release candidate; source review and the institutional code-release decision are needed before all historical scripts and manifests are committed to this public tree.
- On 19 November 2026, a read-only independent GFF3 reannotation audit reproduced **96/96** historical gene/exon/CDS classification labels under both tested deletion-anchor conventions, with **zero discordances**. Reference: *Capra hircus* ARS1.2 RefSeq Annotation Release 102; GFF3 SHA-256 `b702d44f0ebdd638100bc7b93a08163a4e35a61e6651cb398e49d542d5a0b826`.

## Sample provenance and accessions

- Investigator-confirmed matched animal: D01 = 803M = Illumina analysis label 803m.
- NCBI BioProject: [PRJNA1502507](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1502507).
- Publicly resolved BioSample for 803M: [SAMN61979500](https://www.ncbi.nlm.nih.gov/biosample/SAMN61979500) (male; blood; Kalyani, West Bengal; sampled 2026-03-30).
- HiFi SRA accessions cited in the discovery manuscript: SRR41146742 (803M) and SRR41146743 (D863F). Public release status requires confirmation.
- **The 22 Illumina SRR accessions are not verified**. Do not label the HiFi runs as Illumina data. Public linkage searches returned no SRA run IDs at the time of audit. The rest of the cohort's sample-to-library crosswalk and release approval remain pending.

## Publications and releases

- **Historical software DOI (V0.2.0):** https://doi.org/10.5281/zenodo.24368553 — not a genotype-dataset DOI and not a citation for later corrections.
- **Discovery code:** https://github.com/sshil1/bbgoat-sv-scripts
- **22-animal genotype dataset DOI:** pending approved Zenodo deposit.
- Current manuscript: provenance-audited working draft, **not yet authorized for journal submission**. Sampling ethics, authorship, individual Illumina run accessions and institutional data-sharing approval must be documented.

For traceable limitations and release conditions, see [docs/PROVENANCE_AND_RELEASE_AUDIT.md](docs/PROVENANCE_AND_RELEASE_AUDIT.md).

## Reproducibility status

The old repository tree contained a ZIP archive rather than the extracted pipeline; consequently, its V0.2.0 GitHub source archive alone does **not** represent the complete Step 35 workflow. The authors have prepared a corrective source release candidate (v0.2.1), but **no V0.2.1 release or DOI should be claimed until the actual source tree is verified and a new release is published**. Do not rewrite the historical V0.2.0 tag.
