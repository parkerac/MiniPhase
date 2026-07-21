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
  --father-vcf father.vcf.gz \
  --mother-vcf mother.vcf.gz \
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

It can also contain sample-specific input columns. The `bam` column can contain BAM or CRAM paths.

```text
chrom	pos1	ref1	alt1	pos2	ref2	alt2	bam	vcf	sample	father_vcf	mother_vcf	father_sample	mother_sample
chr1	100000	A	G	100120	C	T	sample1.bam	sample1.vcf.gz	SAMPLE1	father1.vcf.gz	mother1.vcf.gz	FATHER1	MOTHER1
chr1	200000	G	A	200180	T	C	sample2.cram	sample2.vcf.gz	SAMPLE2	father2.vcf.gz	mother2.vcf.gz	FATHER2	MOTHER2
chr1	300000	C	A	300150	G	T	sample3.bam	sample3.vcf.gz	SAMPLE3	NA	NA	NA	NA
```

`bam`, `vcf`, `sample`, `father_vcf`, `mother_vcf`, `father_sample`, and `mother_sample` columns override the matching command-line values for that row. This means you can provide common defaults globally and only include columns that vary by sample. `--reference` is shared across all rows and is required when any alignment file is CRAM.

Parent columns can use `NA`, `N/A`, `.`, `None`, `null`, or an empty value when parent data is unavailable for a sample. Those rows skip trio first-pass phasing and use read-backed phasing.

## Trio VCF Phasing

Trio evidence is VCF-only and optional. Provide both parent VCFs to enable it:

```bash
python3 phasing/scripts/phase_nearby_variants.py \
  --reference GRCh38.fa \
  --pairs-tsv variant_pairs.tsv \
  --father-vcf father.vcf.gz \
  --mother-vcf mother.vcf.gz \
  --threads 8 \
  --out phase_results.tsv
```

The script labels a proband ALT allele as `maternal` when the mother carries the ALT and the father is homozygous reference, and as `paternal` for the reverse.

By default, target-variant trio phasing is the first pass. If both target variants have clear parental origins, the script reports `phase` from trio evidence and does not open the proband BAM/CRAM. If trio evidence is not informative, it falls back to read-backed phasing and can still use trio-labeled bridge variants in the local graph.

`bam` is required only when trio first-pass phasing is uninformative or when `--always-run-reads` is used.

Trio-specific output columns include `trio_phase`, `trio_score`, `variant1_origin`, `variant2_origin`, `trio_informative_variants`, `trio_conflicts`, and `trio_status`. The `method` column reports `trio_first_pass` when the proband alignment was skipped and `read_backed` when read-backed phasing was run.

`trio_status` values:

```text
available
missing_parent_vcf
missing_parent_file
disabled
```

Useful trio options:

```text
--trio-weight 10
--min-parent-gq 20
--min-parent-dp 8
--always-run-reads
```

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
