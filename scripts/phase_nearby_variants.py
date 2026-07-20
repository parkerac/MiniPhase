#!/usr/bin/env python3
"""Phase two nearby variants from local BAM/CRAM read evidence without downsampling."""

import argparse
import csv
import json
import math
import multiprocessing as mp
from collections import defaultdict
from dataclasses import dataclass

pysam = None


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    name: str
    target: bool = False


def parse_variant(value, name, target=True):
    parts = value.replace(",", ":").split(":")
    if len(parts) != 4:
        raise ValueError(f"{name} must look like chrom:pos:ref:alt")
    return Variant(parts[0], int(parts[1]), parts[2].upper(), parts[3].upper(), name, target)


def parse_gt(gt):
    if gt is None:
        return None
    alleles = [a for a in gt if a is not None]
    if len(alleles) != 2 or alleles[0] == alleles[1]:
        return None
    return set(alleles)


def same_variant(a, b):
    return a.chrom == b.chrom and a.pos == b.pos and a.ref == b.ref and a.alt == b.alt


def load_bridge_variants(vcf_path, sample, chrom, start, end, targets, max_bridges):
    if not vcf_path:
        return []
    variants = []
    with pysam.VariantFile(vcf_path) as vcf:
        sample = sample or next(iter(vcf.header.samples), None)
        if not sample:
            raise ValueError("VCF has no samples; provide targets only or a genotyped VCF")
        for rec in vcf.fetch(chrom, max(0, start - 1), end):
            if len(rec.ref) > 50 or len(rec.alts or []) != 1:
                continue
            gt = parse_gt(rec.samples[sample].get("GT"))
            if gt != {0, 1}:
                continue
            alt = rec.alts[0]
            variant = Variant(rec.chrom, rec.pos, rec.ref.upper(), alt.upper(), rec.id or f"{rec.chrom}:{rec.pos}:{rec.ref}>{alt}")
            if any(same_variant(variant, target) for target in targets):
                continue
            variants.append(variant)
            if len(variants) >= max_bridges:
                break
    return variants


def insertion_after(aligned_pairs, pair_index, read):
    bases = []
    for qpos, rpos in aligned_pairs[pair_index + 1 :]:
        if rpos is not None:
            break
        if qpos is not None:
            bases.append(read.query_sequence[qpos].upper())
    return "".join(bases)


def call_variant(read, variant, min_baseq):
    if read.is_unmapped or read.query_sequence is None:
        return None
    start = variant.pos - 1
    end = start + len(variant.ref)
    pairs = read.get_aligned_pairs(matches_only=False)
    ref_to_qpos = {rpos: qpos for qpos, rpos in pairs if rpos is not None}
    qpos = [ref_to_qpos.get(pos) for pos in range(start, end)]
    quals = read.query_qualities or []
    observed = "".join("-" if q is None else read.query_sequence[q].upper() for q in qpos)
    baseq = min([quals[q] for q in qpos if q is not None] or [0])

    if len(variant.ref) == len(variant.alt):
        if None in qpos or baseq < min_baseq:
            return None
        if observed == variant.ref:
            return 0, baseq
        if observed == variant.alt:
            return 1, baseq
        return None

    if len(variant.ref) > len(variant.alt) and variant.ref.startswith(variant.alt):
        padded_alt = variant.alt + "-" * (len(variant.ref) - len(variant.alt))
        if baseq >= min_baseq and observed == variant.ref:
            return 0, baseq
        if baseq >= min_baseq and observed == padded_alt:
            return 1, baseq
        return None

    if len(variant.alt) > len(variant.ref) and variant.alt.startswith(variant.ref):
        if None in qpos or baseq < min_baseq or observed != variant.ref:
            return None
        last_ref = end - 1
        pair_index = next((i for i, (q, r) in enumerate(pairs) if r == last_ref and q == qpos[-1]), None)
        inserted = insertion_after(pairs, pair_index, read) if pair_index is not None else ""
        if inserted == "":
            return 0, baseq
        if variant.ref + inserted == variant.alt:
            return 1, baseq
    return None


def add_call(fragment, variant_index, allele, quality):
    old = fragment.get(variant_index)
    if old is None or old[0] == allele and quality > old[1]:
        fragment[variant_index] = (allele, quality)
    elif old[0] != allele:
        fragment[variant_index] = None


def is_cram(path):
    return path.lower().endswith(".cram")


def read_fragments(bam_path, reference, variants, start, end, min_mapq, min_baseq, include_duplicates):
    fragments = defaultdict(dict)
    if is_cram(bam_path) and not reference:
        raise ValueError(f"{bam_path} is a CRAM; provide the shared reference with --reference")
    bam = pysam.AlignmentFile(bam_path, "rc" if is_cram(bam_path) else "rb", reference_filename=reference)
    for read in bam.fetch(variants[0].chrom, max(0, start - 1), end):
        if read.mapping_quality < min_mapq or read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        if read.is_duplicate and not include_duplicates:
            continue
        key = read.query_name
        for i, variant in enumerate(variants):
            if read.reference_start > variant.pos - 1 or read.reference_end is None or read.reference_end < variant.pos:
                continue
            call = call_variant(read, variant, min_baseq)
            if call:
                add_call(fragments[key], i, call[0], min(read.mapping_quality, call[1]))
    bam.close()
    return {name: {i: c for i, c in calls.items() if c is not None} for name, calls in fragments.items()}


def edge_votes(fragments):
    votes = defaultdict(lambda: [0.0, 0.0, 0])
    informative = 0
    for calls in fragments.values():
        items = sorted(calls.items())
        if len(items) < 2:
            continue
        informative += 1
        for a in range(len(items)):
            i, (allele_i, qual_i) = items[a]
            for j, (allele_j, qual_j) in items[a + 1 :]:
                relation = allele_i ^ allele_j
                weight = max(1.0, min(60, qual_i, qual_j) / 10.0)
                votes[(i, j)][relation] += weight
                votes[(i, j)][2] += 1
    return votes, informative


def best_path(votes, source, target, skip_edge=None):
    graph = defaultdict(list)
    for (i, j), (cis, trans, _) in votes.items():
        if skip_edge and (i, j) == skip_edge:
            continue
        if cis == trans:
            continue
        relation = 0 if cis > trans else 1
        support = abs(cis - trans)
        graph[i].append((j, relation, support))
        graph[j].append((i, relation, support))
    best = {(source, 0): (math.inf, [source])}
    queue = [(source, 0)]
    while queue:
        node, parity = queue.pop(0)
        score, path = best[(node, parity)]
        for nxt, edge_parity, support in graph[node]:
            new_parity = parity ^ edge_parity
            new_score = min(score, support)
            state = (nxt, new_parity)
            if new_score > best.get(state, (-1, []))[0]:
                best[state] = (new_score, path + [nxt])
                queue.append(state)
    return best.get((target, 0), (0.0, [])), best.get((target, 1), (0.0, []))


def classify_pair(variants, fragments):
    votes, informative = edge_votes(fragments)
    direct = votes.get((0, 1), [0.0, 0.0, 0])
    cis_path, trans_path = best_path(votes, 0, 1, skip_edge=(0, 1))
    direct_delta = direct[0] - direct[1]
    path_delta = cis_path[0] - trans_path[0]
    score = direct_delta + path_delta
    if abs(score) < 1.0:
        phase = "ambiguous"
    else:
        phase = "cis" if score > 0 else "trans"
    return {
        "phase": phase,
        "score": round(score, 3),
        "direct_cis_weight": round(direct[0], 3),
        "direct_trans_weight": round(direct[1], 3),
        "direct_fragments": direct[2],
        "best_cis_path_score": round(cis_path[0], 3),
        "best_trans_path_score": round(trans_path[0], 3),
        "best_cis_path": ",".join(variants[i].name for i in cis_path[1]),
        "best_trans_path": ",".join(variants[i].name for i in trans_path[1]),
        "informative_fragments": informative,
        "called_fragments": sum(1 for calls in fragments.values() if calls),
    }


def write_fragment_evidence(path, variants, fragments):
    if not path:
        return
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=["fragment", "n_variants", "calls"])
        writer.writeheader()
        for name, calls in sorted(fragments.items()):
            if len(calls) < 2:
                continue
            call_text = ";".join(f"{variants[i].name}:{allele}:q{qual}" for i, (allele, qual) in sorted(calls.items()))
            writer.writerow({"fragment": name, "n_variants": len(calls), "calls": call_text})


def phase_one(args, variant1, variant2):
    if not args.bam:
        raise ValueError("Missing alignment path; provide --bam or a bam column in --pairs-tsv")
    if variant1.chrom != variant2.chrom:
        raise ValueError("Both variants must be on the same chromosome")
    targets = [variant1, variant2]
    span_start = min(variant1.pos, variant2.pos) - args.window
    span_end = max(variant1.pos + len(variant1.ref), variant2.pos + len(variant2.ref)) + args.window
    bridges = load_bridge_variants(args.vcf, args.sample, variant1.chrom, span_start, span_end, targets, args.max_bridges)
    variants = targets + sorted(bridges, key=lambda v: v.pos)
    fragments = read_fragments(args.bam, args.reference, variants, span_start, span_end, args.min_mapq, args.min_baseq, args.include_duplicates)
    result = classify_pair(variants, fragments)
    result.update(
        {
            "variant1": f"{variant1.chrom}:{variant1.pos}:{variant1.ref}:{variant1.alt}",
            "variant2": f"{variant2.chrom}:{variant2.pos}:{variant2.ref}:{variant2.alt}",
            "region": f"{variant1.chrom}:{max(1, span_start)}-{span_end}",
            "bam": args.bam,
            "vcf": args.vcf or "",
            "sample": args.sample or "",
            "bridge_variants": len(bridges),
        }
    )
    return result, variants, fragments


def phase_one_row(item):
    row_index, args, variant1, variant2 = item
    result, _, _ = phase_one(args, variant1, variant2)
    result["row_index"] = row_index
    return result


def row_value(row, key, default=None):
    value = row.get(key)
    return default if value is None or value == "" else value


def pair_rows(args):
    if args.pairs_tsv:
        with open(args.pairs_tsv, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            required = {"chrom", "pos1", "ref1", "alt1", "pos2", "ref2", "alt2"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{args.pairs_tsv} is missing columns: {', '.join(sorted(missing))}")
            for row_index, row in enumerate(reader, start=1):
                row_args = argparse.Namespace(**vars(args))
                row_args.bam = row_value(row, "bam", args.bam)
                row_args.vcf = row_value(row, "vcf", args.vcf)
                row_args.sample = row_value(row, "sample", args.sample)
                yield (
                    row_index,
                    row_args,
                    Variant(row["chrom"], int(row["pos1"]), row["ref1"].upper(), row["alt1"].upper(), "variant1", True),
                    Variant(row["chrom"], int(row["pos2"]), row["ref2"].upper(), row["alt2"].upper(), "variant2", True),
                )
    else:
        yield 1, args, parse_variant(args.variant1, "variant1"), parse_variant(args.variant2, "variant2")


def init_worker():
    global pysam
    import pysam as pysam_module

    pysam = pysam_module


def main():
    global pysam
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", help="Coordinate-sorted, indexed BAM or CRAM; optional in batch mode if pairs TSV has bam")
    parser.add_argument("--reference", help="Shared reference FASTA; required for CRAM and recommended for indels")
    parser.add_argument("--variant1", help="First target as chrom:pos:ref:alt")
    parser.add_argument("--variant2", help="Second target as chrom:pos:ref:alt")
    parser.add_argument("--pairs-tsv", help="Batch TSV with chrom,pos1,ref1,alt1,pos2,ref2,alt2 and optional bam,vcf,sample")
    parser.add_argument("--vcf", help="Optional indexed VCF/BCF of nearby heterozygous bridge variants; can be overridden per batch row")
    parser.add_argument("--sample", help="Sample name in VCF; can be overridden per batch row; defaults to first sample")
    parser.add_argument("--window", type=int, default=1000, help="Bases to fetch around target pair")
    parser.add_argument("--max-bridges", type=int, default=200, help="Maximum nearby heterozygous variants to use")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--include-duplicates", action="store_true")
    parser.add_argument("--threads", type=int, default=1, help="Worker processes for batch mode")
    parser.add_argument("--chunksize", type=int, default=1, help="Pairs assigned per worker task in batch mode")
    parser.add_argument("--out", required=True, help="Output result TSV")
    parser.add_argument("--evidence", help="Optional fragment evidence TSV; only valid for one pair")
    parser.add_argument("--json", help="Optional JSON output; only valid for one pair")
    args = parser.parse_args()

    try:
        import pysam as pysam_module
    except ImportError as exc:
        raise SystemExit("This script requires pysam. Install with: pip install pysam") from exc
    pysam = pysam_module

    if bool(args.pairs_tsv) == bool(args.variant1 or args.variant2):
        raise SystemExit("Provide either --pairs-tsv or both --variant1 and --variant2")
    if not args.pairs_tsv and not (args.variant1 and args.variant2):
        raise SystemExit("Provide both --variant1 and --variant2")
    if not args.pairs_tsv and not args.bam:
        raise SystemExit("Provide --bam for single-pair mode")
    if args.pairs_tsv and (args.evidence or args.json):
        raise SystemExit("--evidence and --json are only supported for single-pair runs")
    if args.threads < 1 or args.chunksize < 1:
        raise SystemExit("--threads and --chunksize must be at least 1")

    rows = []
    last = None
    jobs = list(pair_rows(args))
    if args.pairs_tsv and args.threads > 1:
        with mp.Pool(args.threads, initializer=init_worker) as pool:
            rows = list(pool.imap(phase_one_row, jobs, chunksize=args.chunksize))
    else:
        for row_index, row_args, variant1, variant2 in jobs:
            result, variants, fragments = phase_one(row_args, variant1, variant2)
            result["row_index"] = row_index
            rows.append(result)
            last = result, variants, fragments

    with open(args.out, "w", newline="") as fh:
        fieldnames = list(rows[0].keys()) if rows else []
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    if last and not args.pairs_tsv:
        result, variants, fragments = last
        write_fragment_evidence(args.evidence, variants, fragments)
        if args.json:
            with open(args.json, "w") as fh:
                json.dump(result, fh, indent=2)


if __name__ == "__main__":
    main()
