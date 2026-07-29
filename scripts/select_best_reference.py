#!/usr/bin/env python3
"""Select a reference from full-genome BLAST hits by cumulative bit score."""

import argparse
from collections import defaultdict


def cumulative_bitscores(path):
    scores = defaultdict(float)
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                scores[fields[1]] += float(fields[11])
            except ValueError:
                continue
    return scores


def main():
    parser = argparse.ArgumentParser(
        description="Select the highest cumulative full-genome BLAST bit-score reference."
    )
    parser.add_argument("--blast", required=True, help="BLAST outfmt 6 file with bit score in column 12")
    parser.add_argument("--out", required=True, help="Output file for the selected reference ID")
    parser.add_argument("--scores", required=True, help="Output ranked cumulative score table")
    args = parser.parse_args()

    scores = cumulative_bitscores(args.blast)
    ranked = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    with open(args.scores, "w", encoding="utf-8") as handle:
        handle.write("reference_id\tcumulative_bitscore\n")
        for reference, score in ranked:
            handle.write(f"{reference}\t{score:.6f}\n")

    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write((ranked[0][0] if ranked else "No_Match_Found") + "\n")


if __name__ == "__main__":
    main()
