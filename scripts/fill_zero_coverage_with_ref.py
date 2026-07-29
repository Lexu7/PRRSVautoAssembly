#!/usr/bin/env python3
import argparse
import re
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
    return header, "".join(seq_parts).upper()


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


CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")


def query_to_reference_coordinates(query_length, query_start, target_start, cigar):
    """Map each query coordinate to its homologous reference coordinate."""
    coordinates = [None] * query_length
    query_pos = query_start
    target_pos = target_start
    consumed = 0
    for match in CIGAR_PATTERN.finditer(cigar):
        length = int(match.group(1))
        operation = match.group(2)
        consumed += len(match.group(0))
        if operation in "M=X":
            for _ in range(length):
                if query_pos >= query_length:
                    raise ValueError("CIGAR consumes query coordinates beyond the consensus length")
                coordinates[query_pos] = target_pos
                query_pos += 1
                target_pos += 1
        elif operation in "IS":
            query_pos += length
        elif operation in "DN":
            target_pos += length
        elif operation in "HP":
            continue
        else:
            raise ValueError(f"Unsupported CIGAR operation: {operation}")
    if consumed != len(cigar):
        raise ValueError(f"Invalid CIGAR: {cigar}")
    if query_pos > query_length:
        raise ValueError("CIGAR consumes query coordinates beyond the consensus length")
    return coordinates


def cigar_query_consumed(cigar):
    consumed = 0
    parsed = 0
    for match in CIGAR_PATTERN.finditer(cigar):
        length = int(match.group(1))
        operation = match.group(2)
        parsed += len(match.group(0))
        if operation in "M=XIS":
            consumed += length
    if parsed != len(cigar):
        raise ValueError(f"Invalid CIGAR: {cigar}")
    return consumed


def read_primary_paf(path, query_length):
    primary = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            tags = {field[:5]: field[5:] for field in fields[12:] if len(field) >= 5}
            if tags.get("tp:A:") != "P":
                continue
            cigar = tags.get("cg:Z:")
            if not cigar:
                raise ValueError(f"Primary PAF alignment lacks cg:Z CIGAR: {path}")
            if fields[4] != "+":
                raise ValueError("Reverse-strand reference alignment is not supported for ref filling")
            query_start = int(fields[2])
            query_end = int(fields[3])
            if cigar_query_consumed(cigar) != query_end - query_start:
                raise ValueError(f"PAF query span does not match CIGAR in {path}: {line.rstrip()}")
            primary.append((query_start, query_end, int(fields[7]), cigar))
    if not primary:
        raise ValueError(f"No primary PAF alignment found in {path}")

    coordinates = [None] * query_length
    previous_query_end = 0
    for query_start, query_end, target_start, cigar in sorted(primary):
        if query_start < previous_query_end:
            raise ValueError(f"Overlapping primary PAF query segments in {path}")
        segment_coordinates = query_to_reference_coordinates(
            query_length, query_start, target_start, cigar
        )
        coordinates[query_start:query_end] = segment_coordinates[query_start:query_end]
        previous_query_end = query_end
    return coordinates


def fill_zero_with_reference_coordinates(consensus, best_ref, depths, query_start, target_start, cigar):
    coordinates = query_to_reference_coordinates(len(consensus), query_start, target_start, cigar)
    out = list(consensus)
    filled = 0
    zero_without_reference_coordinate = 0
    nonzero_kept = 0
    for idx, depth in enumerate(depths):
        if depth == 0:
            reference_idx = coordinates[idx]
            if reference_idx is not None and reference_idx < len(best_ref):
                out[idx] = best_ref[reference_idx]
                filled += 1
            else:
                out[idx] = "N"
                zero_without_reference_coordinate += 1
        else:
            nonzero_kept += 1
    return "".join(out), {
        "filled_zero_with_ref_aligned": filled,
        "zero_without_reference_coordinate": zero_without_reference_coordinate,
        "nonzero_kept": nonzero_kept,
    }


def fill_zero_with_ref(consensus, best_ref, depths, alignment_paf):
    coordinates = read_primary_paf(alignment_paf, len(consensus))
    out = list(consensus)
    filled = 0
    zero_without_reference_coordinate = 0
    nonzero_kept = 0
    for idx, depth in enumerate(depths):
        if depth == 0:
            reference_idx = coordinates[idx]
            if reference_idx is not None and reference_idx < len(best_ref):
                out[idx] = best_ref[reference_idx]
                filled += 1
            else:
                out[idx] = "N"
                zero_without_reference_coordinate += 1
        else:
            nonzero_kept += 1
    return "".join(out), filled, zero_without_reference_coordinate, nonzero_kept


def write_fasta(path, sample, seq, line_width=60):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f">{sample}\n")
        for start in range(0, len(seq), line_width):
            handle.write(seq[start : start + line_width] + "\n")


def write_stats(path, original_len, ref_len, output_len, filled, zero_without_ref, nonzero_kept):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("rule\tzero_depth_positions_filled_from_initial_best_ref_by_alignment_coordinate;unmapped_positions_masked_to_N\n")
        handle.write(f"original_length\t{original_len}\n")
        handle.write(f"best_ref_length\t{ref_len}\n")
        handle.write(f"output_length\t{output_len}\n")
        handle.write(f"filled_zero_with_ref\t{filled}\n")
        handle.write(f"zero_without_ref\t{zero_without_ref}\n")
        handle.write(f"filled_zero_with_ref_aligned\t{filled}\n")
        handle.write(f"zero_without_reference_coordinate\t{zero_without_ref}\n")
        handle.write(f"nonzero_kept\t{nonzero_kept}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Create an alternate consensus where zero-depth positions are filled from the initial best reference."
    )
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--best-ref", required=True)
    parser.add_argument("--depth", required=True, help="samtools depth -aa output against the consensus")
    parser.add_argument("--alignment-paf", required=True, help="Primary minimap2 PAF alignment: target=best reference, query=consensus")
    parser.add_argument("--sample", required=True, help="FASTA header for output")
    parser.add_argument("--out", required=True)
    parser.add_argument("--stats", required=True)
    args = parser.parse_args()

    _, consensus = read_single_fasta(args.consensus)
    _, best_ref = read_single_fasta(args.best_ref)
    depths = read_depths(args.depth, len(consensus))
    filled_seq, filled, zero_without_ref, nonzero_kept = fill_zero_with_ref(
        consensus, best_ref, depths, args.alignment_paf
    )

    write_fasta(Path(args.out), args.sample, filled_seq)
    write_stats(
        Path(args.stats),
        len(consensus),
        len(best_ref),
        len(filled_seq),
        filled,
        zero_without_ref,
        nonzero_kept,
    )


if __name__ == "__main__":
    main()
