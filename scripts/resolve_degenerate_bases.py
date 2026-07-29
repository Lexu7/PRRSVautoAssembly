#!/usr/bin/env python3
import argparse
from collections import Counter
from pathlib import Path


DEGENERATE_BASES = set("RYSWKMBDHV")
ACGT = set("ACGT")
IUPAC_ALLOWED = {
    "R": set("AG"),
    "Y": set("CT"),
    "S": set("GC"),
    "W": set("AT"),
    "K": set("GT"),
    "M": set("AC"),
    "B": set("CGT"),
    "D": set("AGT"),
    "H": set("ACT"),
    "V": set("ACG"),
}
BASES_TO_IUPAC = {
    frozenset("A"): "A",
    frozenset("C"): "C",
    frozenset("G"): "G",
    frozenset("T"): "T",
    frozenset("AG"): "R",
    frozenset("CT"): "Y",
    frozenset("GC"): "S",
    frozenset("AT"): "W",
    frozenset("GT"): "K",
    frozenset("AC"): "M",
    frozenset("CGT"): "B",
    frozenset("AGT"): "D",
    frozenset("ACT"): "H",
    frozenset("ACG"): "V",
}


def read_single_fasta(path):
    header = None
    seq_parts = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    raise ValueError(f"Expected one FASTA record, found multiple in {path}")
                header = line[1:].strip()
            else:
                seq_parts.append(line.strip())
    if header is None:
        raise ValueError(f"No FASTA record found in {path}")
    return header, "".join(seq_parts)


def write_fasta(path, header, seq, line_width=60):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f">{header}\n")
        for start in range(0, len(seq), line_width):
            handle.write(seq[start : start + line_width] + "\n")


def skip_numbered_payload(bases, index):
    index += 1
    start = index
    while index < len(bases) and bases[index].isdigit():
        index += 1
    if start == index:
        return index
    length = int(bases[start:index])
    return index + length


def count_mpileup_bases(bases, ref_base):
    counts = Counter()
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
            i = skip_numbered_payload(bases, i)
            continue
        if char in ".,":  # Match to the reference base.
            if ref_base in ACGT:
                counts[ref_base] += 1
            i += 1
            continue
        base = char.upper()
        if base in ACGT:
            counts[base] += 1
        i += 1
    return counts


def load_counts_by_position(mpileup_path, seq_len):
    by_pos = {}
    with open(mpileup_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if not 1 <= pos <= seq_len:
                continue
            ref_base = parts[2]
            bases = parts[4]
            by_pos[pos] = count_mpileup_bases(bases, ref_base)
    return by_pos


def top_bases(counts):
    if not counts:
        return [], 0
    top_count = max(counts.values())
    return sorted(base for base, count in counts.items() if count == top_count), top_count


def resolve_degenerate(seq, counts_by_pos):
    resolved = list(seq)
    total = 0
    to_acgt = 0
    to_iupac = 0
    to_n = 0
    no_pileup = 0
    incompatible_top = 0
    tied_top = 0
    unrepresentable_tie = 0

    for idx, base in enumerate(seq):
        base = base.upper()
        if base not in DEGENERATE_BASES:
            continue
        total += 1
        counts = counts_by_pos.get(idx + 1, Counter())
        depth = sum(counts.values())
        if depth == 0:
            no_pileup += 1
            resolved[idx] = "N"
            to_n += 1
            continue
        top, _ = top_bases(counts)
        if len(top) > 1:
            tied_top += 1
            iupac = BASES_TO_IUPAC.get(frozenset(top))
            if iupac:
                resolved[idx] = iupac
                to_iupac += 1
            else:
                # Four-way ties can only be represented by N. Keep a concrete
                # observed base instead and report the tie.
                resolved[idx] = top[0]
                to_acgt += 1
                unrepresentable_tie += 1
            continue
        top_base = top[0]
        if top_base not in IUPAC_ALLOWED[base]:
            incompatible_top += 1
        resolved[idx] = top_base
        to_acgt += 1

    stats = {
        "degenerate_total": total,
        "resolved_to_acgt": to_acgt,
        "resolved_to_iupac": to_iupac,
        "masked_to_n": to_n,
        "no_pileup": no_pileup,
        "incompatible_top": incompatible_top,
        "tied_top": tied_top,
        "unrepresentable_tie": unrepresentable_tie,
    }
    return "".join(resolved), stats


def write_stats(path, stats):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("rule\tuse_observed_top_base_only_zero_depth_to_N\n")
        for key in (
            "degenerate_total",
            "resolved_to_acgt",
            "resolved_to_iupac",
            "masked_to_n",
            "no_pileup",
            "incompatible_top",
            "tied_top",
            "unrepresentable_tie",
        ):
            handle.write(f"{key}\t{stats[key]}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Resolve IUPAC degenerate bases in a FASTA using samtools mpileup support."
    )
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--mpileup", default="")
    parser.add_argument("--no-degenerate", action="store_true", help="Write unchanged FASTA and zero statistics without reading mpileup.")
    parser.add_argument("--min-dp", type=int, default=None, help="Ignored; only zero A/C/G/T observations become N.")
    parser.add_argument("--min-frac", type=float, default=None, help="Ignored; kept for backward compatibility.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--stats", required=True)
    args = parser.parse_args()

    header, seq = read_single_fasta(args.fasta)
    if args.no_degenerate:
        stats = {
            "degenerate_total": 0,
            "resolved_to_acgt": 0,
            "resolved_to_iupac": 0,
            "masked_to_n": 0,
            "no_pileup": 0,
            "incompatible_top": 0,
            "tied_top": 0,
            "unrepresentable_tie": 0,
        }
        write_fasta(Path(args.out), header, seq)
        write_stats(Path(args.stats), stats)
        return
    if not args.mpileup:
        parser.error("--mpileup is required unless --no-degenerate is used")
    counts_by_pos = load_counts_by_position(args.mpileup, len(seq))
    resolved, stats = resolve_degenerate(seq, counts_by_pos)
    write_fasta(Path(args.out), header, resolved)
    write_stats(Path(args.stats), stats)


if __name__ == "__main__":
    main()
