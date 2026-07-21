#!/usr/bin/env python3
"""Phase two nearby variants from local BAM/CRAM read evidence without downsampling."""

import argparse
import csv
import json
import math
import multiprocessing as mp
import os
from collections import defaultdict
from dataclasses import dataclass

pysam = None
MISSING_VALUES = {"", ".", "NA", "N/A", "NONE", "NULL", "NAN"}


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


def alt_status(sample_call, alt_index, min_gq, min_dp):
    gt = sample_call.get("GT") if sample_call else None
    if not gt or any(a is None for a in gt):
        return "missing"
    if "GQ" in sample_call and sample_call["GQ"] is not None and sample_call["GQ"] < min_gq:
        return "low_quality"
    if "DP" in sample_call and sample_call["DP"] is not None and sample_call["DP"] < min_dp:
        return "low_quality"
    alleles = set(gt)
    if alt_index in alleles:
        return "has_alt"
    if alleles == {0}:
        return "no_alt"
    return "other_alt"


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


def load_parent_statuses(vcf_path, sample, chrom, start, end, variants, min_gq, min_dp):
    statuses = {i: "missing" for i in range(len(variants))}
    if not vcf_path:
        return statuses, sample or ""
    wanted = {(v.chrom, v.pos, v.ref, v.alt): i for i, v in enumerate(variants)}
    with pysam.VariantFile(vcf_path) as vcf:
        sample = sample or next(iter(vcf.header.samples), None)
        if not sample:
            return statuses, ""
        if sample not in vcf.header.samples:
            raise ValueError(f"{sample} is not present in {vcf_path}")
        for rec in vcf.fetch(chrom, max(0, start - 1), end):
            if not rec.alts:
                continue
            for alt_index, alt in enumerate(rec.alts, start=1):
                key = (rec.chrom, rec.pos, rec.ref.upper(), alt.upper())
                if key in wanted:
                    statuses[wanted[key]] = alt_status(rec.samples[sample], alt_index, min_gq, min_dp)
    return statuses, sample


def infer_origin(mother_status, father_status):
    if mother_status == "low_quality" or father_status == "low_quality":
        return "low_quality"
    if mother_status == "has_alt" and father_status == "no_alt":
        return "maternal"
    if father_status == "has_alt" and mother_status == "no_alt":
        return "paternal"
    if mother_status == "no_alt" and father_status == "no_alt":
        return "conflicting"
    return "uninformative"


def trio_phase_from_origins(origin1, origin2):
    if origin1 == "conflicting" or origin2 == "conflicting":
        return "conflicting"
    if origin1 not in {"maternal", "paternal"} or origin2 not in {"maternal", "paternal"}:
        return "ambiguous"
    return "cis" if origin1 == origin2 else "trans"


def trio_status(args):
    if args.trio_weight <= 0:
        return "disabled"
    if not args.mother_vcf or not args.father_vcf:
        return "missing_parent_vcf"
    if not os.path.exists(args.mother_vcf) or not os.path.exists(args.father_vcf):
        return "missing_parent_file"
    return "available"


def trio_edges(args, variants, chrom, start, end):
    votes = defaultdict(lambda: [0.0, 0.0, 0])
    origins = {i: "not_tested" for i in range(len(variants))}
    if trio_status(args) != "available":
        return votes, origins, 0, 0
    mother, mother_sample = load_parent_statuses(args.mother_vcf, args.mother_sample, chrom, start, end, variants, args.min_parent_gq, args.min_parent_dp)
    father, father_sample = load_parent_statuses(args.father_vcf, args.father_sample, chrom, start, end, variants, args.min_parent_gq, args.min_parent_dp)
    args.mother_sample = mother_sample
    args.father_sample = father_sample
    for i in range(len(variants)):
        origins[i] = infer_origin(mother[i], father[i])
    known = [(i, origin) for i, origin in origins.items() if origin in {"maternal", "paternal"}]
    for a in range(len(known)):
        i, origin_i = known[a]
        for j, origin_j in known[a + 1 :]:
            relation = 0 if origin_i == origin_j else 1
            votes[(i, j)][relation] += args.trio_weight
            votes[(i, j)][2] += 1
    target_origins = [origins.get(0, "not_tested"), origins.get(1, "not_tested")]
    conflicts = sum(1 for origin in target_origins if origin == "conflicting")
    return votes, origins, len(known), conflicts


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


def merge_votes(*vote_sets):
    merged = defaultdict(lambda: [0.0, 0.0, 0])
    for votes in vote_sets:
        for edge, values in votes.items():
            merged[edge][0] += values[0]
            merged[edge][1] += values[1]
            merged[edge][2] += values[2]
    return merged


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


def classify_pair(variants, fragments, trio_votes=None):
    read_votes, informative = edge_votes(fragments)
    trio_votes = trio_votes or {}
    votes = merge_votes(read_votes, trio_votes)
    direct_read = read_votes.get((0, 1), [0.0, 0.0, 0])
    direct_trio = trio_votes.get((0, 1), [0.0, 0.0, 0])
    direct = votes.get((0, 1), [0.0, 0.0, 0])
    read_cis_path, read_trans_path = best_path(read_votes, 0, 1, skip_edge=(0, 1))
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
        "direct_cis_weight": round(direct_read[0], 3),
        "direct_trans_weight": round(direct_read[1], 3),
        "direct_fragments": direct_read[2],
        "trio_direct_cis_weight": round(direct_trio[0], 3),
        "trio_direct_trans_weight": round(direct_trio[1], 3),
        "trio_score": round(direct_trio[0] - direct_trio[1], 3),
        "best_read_cis_path_score": round(read_cis_path[0], 3),
        "best_read_trans_path_score": round(read_trans_path[0], 3),
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


def add_metadata(result, args, variant1, variant2, span_start, span_end, origins, trio_informative, trio_conflicts, bridge_count, method, status):
    variant1_origin = origins.get(0, "not_tested")
    variant2_origin = origins.get(1, "not_tested")
    result.update(
        {
            "variant1": f"{variant1.chrom}:{variant1.pos}:{variant1.ref}:{variant1.alt}",
            "variant2": f"{variant2.chrom}:{variant2.pos}:{variant2.ref}:{variant2.alt}",
            "region": f"{variant1.chrom}:{max(1, span_start)}-{span_end}",
            "bam": args.bam,
            "vcf": args.vcf or "",
            "sample": args.sample or "",
            "father_vcf": args.father_vcf or "",
            "mother_vcf": args.mother_vcf or "",
            "father_sample": args.father_sample or "",
            "mother_sample": args.mother_sample or "",
            "variant1_origin": variant1_origin,
            "variant2_origin": variant2_origin,
            "trio_phase": trio_phase_from_origins(variant1_origin, variant2_origin),
            "trio_informative_variants": trio_informative,
            "trio_conflicts": trio_conflicts,
            "trio_status": status,
            "bridge_variants": bridge_count,
            "method": method,
        }
    )
    return result


def phase_one(args, variant1, variant2):
    if variant1.chrom != variant2.chrom:
        raise ValueError("Both variants must be on the same chromosome")
    targets = [variant1, variant2]
    span_start = min(variant1.pos, variant2.pos) - args.window
    span_end = max(variant1.pos + len(variant1.ref), variant2.pos + len(variant2.ref)) + args.window

    status = trio_status(args)
    target_trio_votes, target_origins, target_trio_informative, target_trio_conflicts = trio_edges(args, targets, variant1.chrom, span_start, span_end)
    target_trio_phase = trio_phase_from_origins(target_origins.get(0, "not_tested"), target_origins.get(1, "not_tested"))
    if target_trio_phase in {"cis", "trans"} and not args.always_run_reads:
        result = classify_pair(targets, {}, target_trio_votes)
        result["phase"] = target_trio_phase
        return add_metadata(result, args, variant1, variant2, span_start, span_end, target_origins, target_trio_informative, target_trio_conflicts, 0, "trio_first_pass", status), targets, {}

    if not args.bam:
        raise ValueError("Trio first pass was not informative; provide --bam or a bam column in --pairs-tsv for read-backed phasing")
    bridges = load_bridge_variants(args.vcf, args.sample, variant1.chrom, span_start, span_end, targets, args.max_bridges)
    variants = targets + sorted(bridges, key=lambda v: v.pos)
    trio_vote_set, origins, trio_informative, trio_conflicts = trio_edges(args, variants, variant1.chrom, span_start, span_end)
    fragments = read_fragments(args.bam, args.reference, variants, span_start, span_end, args.min_mapq, args.min_baseq, args.include_duplicates)
    result = classify_pair(variants, fragments, trio_vote_set)
    result = add_metadata(result, args, variant1, variant2, span_start, span_end, origins, trio_informative, trio_conflicts, len(bridges), "read_backed", status)
    return result, variants, fragments


def phase_one_row(item):
    row_index, args, variant1, variant2 = item
    result, _, _ = phase_one(args, variant1, variant2)
    result["row_index"] = row_index
    return result


def row_value(row, key, default=None, missing_overrides=False):
    if key not in row:
        return default
    value = row.get(key)
    if value is None or value.strip().upper() in MISSING_VALUES:
        return None if missing_overrides else default
    return value


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
                row_args.father_vcf = row_value(row, "father_vcf", args.father_vcf, missing_overrides=True)
                row_args.mother_vcf = row_value(row, "mother_vcf", args.mother_vcf, missing_overrides=True)
                row_args.father_sample = row_value(row, "father_sample", args.father_sample, missing_overrides=True)
                row_args.mother_sample = row_value(row, "mother_sample", args.mother_sample, missing_overrides=True)
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
    parser.add_argument("--bam", help="Coordinate-sorted, indexed BAM or CRAM; required when trio first pass is uninformative")
    parser.add_argument("--reference", help="Shared reference FASTA; required for CRAM and recommended for indels")
    parser.add_argument("--variant1", help="First target as chrom:pos:ref:alt")
    parser.add_argument("--variant2", help="Second target as chrom:pos:ref:alt")
    parser.add_argument("--pairs-tsv", help="Batch TSV with target variants and optional sample/proband/parent input columns")
    parser.add_argument("--vcf", help="Optional indexed VCF/BCF of nearby heterozygous bridge variants; can be overridden per batch row")
    parser.add_argument("--sample", help="Sample name in VCF; can be overridden per batch row; defaults to first sample")
    parser.add_argument("--father-vcf", help="Optional father VCF/BCF for trio phasing; can be overridden per batch row")
    parser.add_argument("--mother-vcf", help="Optional mother VCF/BCF for trio phasing; can be overridden per batch row")
    parser.add_argument("--father-sample", help="Father sample name; can be overridden per batch row; defaults to first sample")
    parser.add_argument("--mother-sample", help="Mother sample name; can be overridden per batch row; defaults to first sample")
    parser.add_argument("--trio-weight", type=float, default=10.0, help="Weight added by each informative trio relationship")
    parser.add_argument("--min-parent-gq", type=float, default=20.0, help="Minimum parent GQ when present")
    parser.add_argument("--min-parent-dp", type=float, default=8.0, help="Minimum parent DP when present")
    parser.add_argument("--always-run-reads", action="store_true", help="Run proband read-backed phasing even when target trio phasing is clear")
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
