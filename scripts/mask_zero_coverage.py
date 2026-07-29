#!/usr/bin/env python3
import argparse
from pathlib import Path


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


def read_depths(path, seq_len):
    depths = [0] * seq_len
    seen = 0
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                pos = int(parts[1])
                depth = int(parts[2])
            except ValueError:
                continue
            if 1 <= pos <= seq_len:
                depths[pos - 1] = depth
                seen += 1
    if seq_len and seen == 0:
        raise ValueError(f"No usable depth rows found in {path}")
    return depths


def trim_and_mask_zero_coverage(seq, depths):
    seq_len = len(seq)
    left = 0
    while left < seq_len and depths[left] == 0:
        left += 1

    if left == seq_len:
        return "", left, 0, 0

    right = seq_len - 1
    while right >= left and depths[right] == 0:
        right -= 1

    trimmed_seq = list(seq[left : right + 1])
    trimmed_depths = depths[left : right + 1]
    internal_zero = 0
    for idx, depth in enumerate(trimmed_depths):
        if depth == 0:
            trimmed_seq[idx] = "N"
            internal_zero += 1

    trim_right = seq_len - right - 1
    return "".join(trimmed_seq), left, trim_right, internal_zero


def write_fasta(path, sample, seq, line_width=60):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f">{sample}\n")
        if not seq:
            handle.write("\n")
            return
        for start in range(0, len(seq), line_width):
            handle.write(seq[start : start + line_width] + "\n")


def write_stats(path, original_len, final_len, trim_left, trim_right, internal_zero):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f"original_length\t{original_len}\n")
        handle.write(f"final_length\t{final_len}\n")
        handle.write(f"trim_left\t{trim_left}\n")
        handle.write(f"trim_right\t{trim_right}\n")
        handle.write(f"internal_zero_masked\t{internal_zero}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Trim terminal zero-coverage bases and mask internal zero-coverage bases as N."
    )
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--depth", required=True, help="samtools depth -aa output")
    parser.add_argument("--sample", required=True, help="FASTA header for output")
    parser.add_argument("--out", required=True)
    parser.add_argument("--stats", required=True)
    args = parser.parse_args()

    _, seq = read_single_fasta(args.fasta)
    depths = read_depths(args.depth, len(seq))
    masked_seq, trim_left, trim_right, internal_zero = trim_and_mask_zero_coverage(seq, depths)

    out_path = Path(args.out)
    stats_path = Path(args.stats)
    write_fasta(out_path, args.sample, masked_seq)
    write_stats(stats_path, len(seq), len(masked_seq), trim_left, trim_right, internal_zero)


if __name__ == "__main__":
    main()
