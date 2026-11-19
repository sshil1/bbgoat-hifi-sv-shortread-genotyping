# Production, annotation and release audit

Audit date: 2026-11-19 (Asia/Kolkata). This record describes verified evidence, not laboratory approval or a new release.

## Source-backed computational status

- The 22 × 96 source panel yielded 2,112 unique sample–locus attempts: FT-PASS=1,645; FT-PASS with depth ≥10=1,632; FT20 loci=52; conservative loci=47; genotype-variable conservative loci=39.
- Source-panel selection: 48 deletions and 48 insertions, 12 each per class in each of four size strata; seed 20261117 recorded in selection summary.
- Paragraph environment: Bioconda Paragraph 2.3; Python 3.12.3; pysam 0.22.1. For 21 new samples, 21/21 original commands, manifests and calibration files were captured; their source-panel SHA-256 matches `c3ca2b5c24110ff74029d7134d401c4faff561fe02d5e986278e85b6b8f01ad3`. All 21 logs report 96 graphs completed; 11 read-limit warnings across six runs.
- The 803m pilot reused completed batches of 72 + 6 + 18 graphs; another 4-graph run was smoke testing. Original 803m invocations not recovered. No genotyping rerun is necessary.
- Cross-platform comparison across 96 loci: 50 exact genotype matches, 13 discordant calls, 17 FT non-PASS, 16 otherwise non-evaluable. The evaluable exact-concordance fraction is 50/63 (79.37%), not accuracy against a truth set.
- Source GFF3: NCBI ARS1.2 RefSeq Annotation Release 102; SHA-256 `b702d44f0ebdd638100bc7b93a08163a4e35a61e6651cb398e49d542d5a0b826`.
- Independent read-only Step 36 GFF3 annotation comparison: total loci 96; 96 matched anchor-excluded; 96 matched anchor-included; 0 differing between conventions; 0 discordant under both. The original annotation code was not recovered.

## Identity and accession evidence

The investigator confirms D01 = 803M = 803m. NCBI E-utilities links the male Black Bengal goat 803M to BioSample SAMN61979500, BioProject PRJNA1502507, collection 2026-03-30 in Kalyani, blood. The 22 individual Illumina SRA runs and original 22-sample library-chain-of-custody table have not been recovered; public SRA ELink returned zero records. This is **not** evidence that the datasets have not been submitted privately or embargoed.

## Pending documentary gates

1. Sampling and animal ethics approval/exemption for the present 22-animal original research.
2. Verified animal → extraction → Illumina library → BAM → BioSample/SRA mapping, independent of filenames.
3. Institutional/coauthor approval and data-sharing permissions, especially for unpublished source alleles and individual-level genotype data.
4. Final affiliations, author order, contribution statement and funding/declarations.
5. Release of a new, verified GitHub source tree under an appropriately new tag. Do not alter the Zenodo V0.2.0 snapshot.
