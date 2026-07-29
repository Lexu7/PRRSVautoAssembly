#!/usr/bin/env python3
import argparse
from collections import Counter


HEADER = [
    "sample",
    "chrom",
    "pos",
    "ref",
    "mpileup_depth",
    "acgt_depth",
    "top_base",
    "top_base_count",
    "top_base_frac_acgt",
    "top_base_frac_total",
    "del_placeholder_count",
    "del_placeholder_frac_total",
    "del_event_count",
    "del_event_frac_total",
    "ins_event_count",
    "ins_event_frac_total",
    "top_insertions",
    "top_deletions",
    "flags",
    "base_counts",
]


def parse_bases(bases, ref_base):
    counts = Counter()
    insertions = Counter()
    deletions = Counter()
    i = 0
    ref_base = ref_base.upper()
    while i < len(bases):
        char = bases[i]
        if char == "^":
            i += 2
            continue
        if char == "$":
            i += 1
            continue
        if char in "+-":
            sign = char
            i += 1
            start = i
            while i < len(bases) and bases[i].isdigit():
                i += 1
            length = int(bases[start:i] or 0)
            payload = bases[i : i + length].upper()
            i += length
            if sign == "+":
                insertions[payload] += 1
            else:
                deletions[payload] += 1
            continue
        if char in ".,":
            if ref_base in "ACGT":
                counts[ref_base] += 1
            i += 1
            continue
        if char == "*":
            counts["DEL"] += 1
            i += 1
            continue
        base = char.upper()
        if base in "ACGTN":
            counts[base] += 1
        i += 1
    return counts, insertions, deletions


def top_base_stats(counts):
    acgt = {base: counts.get(base, 0) for base in "ACGT"}
    acgt_depth = sum(acgt.values())
    if acgt_depth == 0:
        return acgt_depth, "NA", 0, 0.0
    top_count = max(acgt.values())
    top = sorted(base for base, count in acgt.items() if count == top_count and count > 0)
    top_base = ",".join(top)
    return acgt_depth, top_base, top_count, top_count / acgt_depth


def format_counts(counts):
    return ",".join(f"{key}:{counts.get(key, 0)}" for key in ["A", "C", "G", "T", "N", "DEL"])


def format_top_indels(counter, sign, limit=3):
    if not counter:
        return "."
    return ",".join(f"{sign}{len(seq)}{seq}:{count}" for seq, count in counter.most_common(limit))


def analyze_line(parts, sample, thresholds):
    chrom, pos, ref, depth_text, bases = parts[:5]
    depth = int(depth_text)
    counts, insertions, deletions = parse_bases(bases, ref)
    acgt_depth, top_base, top_count, top_frac_acgt = top_base_stats(counts)
    top_frac_total = (top_count / depth) if depth else 0.0
    del_placeholder = counts.get("DEL", 0)
    del_placeholder_frac = (del_placeholder / depth) if depth else 0.0
    ins_count = sum(insertions.values())
    del_event_count = sum(deletions.values())
    ins_frac = (ins_count / depth) if depth else 0.0
    del_event_frac = (del_event_count / depth) if depth else 0.0

    flags = []
    if depth < thresholds["min_depth"]:
        flags.append("LOW_DEPTH")
    if del_placeholder_frac >= thresholds["del_frac"]:
        flags.append("DEL_DOMINANT")
    if acgt_depth and top_frac_acgt < thresholds["base_frac"]:
        flags.append("MIXED_BASE")
    if ins_frac >= thresholds["ins_frac"]:
        flags.append("INS_SUPPORTED")
    if del_event_count and del_event_frac >= thresholds["del_event_frac"]:
        flags.append("DEL_EVENT_SUPPORTED")

    if not flags:
        return None

    return {
        "sample": sample,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "mpileup_depth": str(depth),
        "acgt_depth": str(acgt_depth),
        "top_base": top_base,
        "top_base_count": str(top_count),
        "top_base_frac_acgt": f"{top_frac_acgt:.4f}",
        "top_base_frac_total": f"{top_frac_total:.4f}",
        "del_placeholder_count": str(del_placeholder),
        "del_placeholder_frac_total": f"{del_placeholder_frac:.4f}",
        "del_event_count": str(del_event_count),
        "del_event_frac_total": f"{del_event_frac:.4f}",
        "ins_event_count": str(ins_count),
        "ins_event_frac_total": f"{ins_frac:.4f}",
        "top_insertions": format_top_indels(insertions, "+"),
        "top_deletions": format_top_indels(deletions, "-"),
        "flags": "|".join(flags),
        "base_counts": format_counts(counts),
    }


def write_summary(path, sample, rows):
    flag_counts = Counter()
    for row in rows:
        for flag in row["flags"].split("|"):
            flag_counts[flag] += 1
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f"sample\t{sample}\n")
        handle.write(f"flagged_sites\t{len(rows)}\n")
        for key in ["LOW_DEPTH", "DEL_DOMINANT", "MIXED_BASE", "INS_SUPPORTED", "DEL_EVENT_SUPPORTED"]:
            handle.write(f"{key.lower()}\t{flag_counts.get(key, 0)}\n")


def main():
    parser = argparse.ArgumentParser(description="Flag indel-complex and mixed-base sites from samtools mpileup.")
    parser.add_argument("--mpileup", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--min-depth", type=int, default=20)
    parser.add_argument("--base-frac", type=float, default=0.80)
    parser.add_argument("--del-frac", type=float, default=0.30)
    parser.add_argument("--ins-frac", type=float, default=0.20)
    parser.add_argument("--del-event-frac", type=float, default=0.20)
    args = parser.parse_args()

    thresholds = {
        "min_depth": args.min_depth,
        "base_frac": args.base_frac,
        "del_frac": args.del_frac,
        "ins_frac": args.ins_frac,
        "del_event_frac": args.del_event_frac,
    }
    rows = []
    with open(args.mpileup, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            row = analyze_line(parts, args.sample, thresholds)
            if row:
                rows.append(row)

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("\t".join(HEADER) + "\n")
        for row in rows:
            handle.write("\t".join(row[column] for column in HEADER) + "\n")
    write_summary(args.summary, args.sample, rows)


if __name__ == "__main__":
    main()
