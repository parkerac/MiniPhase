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
chrom	pos1	ref1	alt1	pos2	ref2	alt2	bam	vcf	sample
chr1	100000	A	G	100120	C	T	sample1.bam	sample1.vcf.gz	SAMPLE1
chr1	200000	G	A	200180	T	C	sample2.cram	sample2.vcf.gz	SAMPLE2
```

`bam`, `vcf`, and `sample` columns override the matching command-line values for that row. This means you can provide common defaults globally and only include columns that vary by sample. `bam` is required either as `--bam` or as a `bam` column. `--reference` is shared across all rows and is required when any alignment file is CRAM.

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
