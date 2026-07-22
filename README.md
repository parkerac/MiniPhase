# Nearby Variant Phaser

`phase_nearby_variants.py` phases two nearby variants from an indexed BAM/CRAM without read downsampling.

It uses:

- reads that directly overlap both target variants
- paired-end fragments, merged by read name
- optional nearby heterozygous bridge variants from an indexed VCF/BCF

## Install dependency

```bash
pip install pysam
```

## Single Pair

```bash
python3 phasing/scripts/phase_nearby_variants.py \
  --bam sample.bam \
  --reference GRCh38.fa \
  --variant1 chr1:100000:A:G \
  --variant2 chr1:100120:C:T \
  --vcf sample.vcf.gz \
  --sample SAMPLE_ID \
  --father-bam father.bam \
  --mother-bam mother.bam \
  --out phase_result.tsv \
  --evidence phase_fragments.tsv \
  --json phase_result.json
```

`phase` is reported as:

- `cis`: alternate alleles are inferred on the same haplotype
- `trans`: alternate alleles are inferred on opposite haplotypes
- `ambiguous`: evidence is absent or balanced

## Batch Mode

```bash
python3 phasing/scripts/phase_nearby_variants.py \
  --reference GRCh38.fa \
  --pairs-tsv variant_pairs.tsv \
  --threads 8 \
  --out phase_results.tsv
```

The batch TSV must contain:

```text
chrom	pos1	ref1	alt1	pos2	ref2	alt2
```

It can also contain sample-specific input columns. The `bam`, `father_bam`, and `mother_bam` columns can contain BAM or CRAM paths.

```text
chrom	pos1	ref1	alt1	pos2	ref2	alt2	bam	reference	vcf	sample	father_bam	mother_bam
chr1	100000	A	G	100120	C	T	sample1.bam	GRCh38.fa	sample1.vcf.gz	SAMPLE1	father1.bam	mother1.bam
chr1	200000	G	A	200180	T	C	sample2.cram	sample2.fa	sample2.vcf.gz	SAMPLE2	father2.cram	mother2.cram
chr1	300000	C	A	300150	G	T	sample3.bam	GRCh38.fa	sample3.vcf.gz	SAMPLE3	NA	NA
```

`bam`, `reference`, `vcf`, `sample`, `father_bam`, and `mother_bam` columns override the matching command-line values for that row. This means you can provide common defaults globally and only include columns that vary by sample. A reference FASTA is required for CRAM inputs and recommended for indels.

Reference paths are canonicalized internally. In batch mode, rows are grouped by reference for processing and then written back in the original input order. Each worker process keeps one open FASTA handle per reference path, which helps avoid repeated reference setup when many rows share the same FASTA.

Contig names are corrected at fetch time for common reference/header mismatches. For example, a target listed as `chr10` can still be fetched from a BAM/CRAM, VCF, or FASTA-indexed setup that uses `10`, and mitochondrial aliases `M`, `MT`, `chrM`, and `chrMT` are handled similarly.

Parent BAM columns can use `NA`, `N/A`, `.`, `None`, `null`, or an empty value when parent data is unavailable for a sample. Those rows still run proband read-backed phasing.

## Phasing Order

The script first tries proband read-backed phasing using direct overlapping reads, paired-end fragments, and optional nearby heterozygous bridge variants from the proband VCF.

If the proband read-backed result is `ambiguous`, the script optionally checks parent BAM/CRAM files at the two target variants only:

```bash
python3 phasing/scripts/phase_nearby_variants.py \
  --reference GRCh38.fa \
  --pairs-tsv variant_pairs.tsv \
  --father-bam father.bam \
  --mother-bam mother.bam \
  --threads 8 \
  --out phase_results.tsv
```

Parent BAM rescue labels a proband ALT allele as `maternal` when the mother has ALT-supporting reads and the father is consistent with reference, and as `paternal` for the reverse. If both target origins are clear, the final `phase` is rescued as `cis` or `trans`.

Useful parent BAM rescue options:

```text
--min-parent-bam-depth 8
--min-parent-bam-alt-depth 3
--min-parent-bam-alt-frac 0.2
--max-parent-bam-alt-depth-for-ref 1
--max-parent-bam-alt-frac-for-ref 0.05
```

The `method` column reports `read_backed`, `parent_bam_rescue`, or `read_backed_parent_bam_uninformative`.

Parent-BAM-specific output columns include `parent_bam_phase`, `variant1_mother_bam_status`, `variant1_father_bam_status`, `variant2_mother_bam_status`, `variant2_father_bam_status`, and matching REF/ALT depth and ALT fraction columns. Parent BAM statuses are `has_alt`, `no_alt`, `low_depth`, `ambiguous`, `missing`, `missing_file`, or `missing_reference`.

`--threads` controls the number of worker processes. Output row order matches the input TSV order.

For very large batches, increasing `--chunksize` can reduce multiprocessing overhead:

```bash
python3 phasing/scripts/phase_nearby_variants.py \
  --reference GRCh38.fa \
  --pairs-tsv variant_pairs.tsv \
  --threads 8 \
  --chunksize 20 \
  --out phase_results.tsv
```
