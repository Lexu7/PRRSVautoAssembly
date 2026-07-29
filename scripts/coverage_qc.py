#!/usr/bin/env python3
import argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", required=True, help="samtools depth -aa output: chrom pos depth")
    ap.add_argument("--low-dp", type=int, default=10, help="depth < low_dp considered low coverage")
    ap.add_argument("--low-frac", type=float, default=0.30, help="warn if low coverage fraction >= this")
    ap.add_argument("--zero-frac", type=float, default=0.05, help="warn if zero coverage fraction >= this")
    ap.add_argument("-o", "--out", required=True, help="output key\\tvalue file")
    args = ap.parse_args()

    total = 0
    zero = 0
    low = 0
    depth_sum = 0

    with open(args.depth) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            d = int(parts[2])
            total += 1
            depth_sum += d
            if d == 0:
                zero += 1
            if d < args.low_dp:
                low += 1

    avg_depth = (depth_sum / total) if total else 0.0
    low_frac = (low / total) if total else 0.0
    zero_frac = (zero / total) if total else 0.0

    warns = []
    if low_frac >= args.low_frac:
        warns.append("LOW_COVERAGE")
    if zero_frac >= args.zero_frac:
        warns.append("ZERO_COVERAGE")

    status = "OK" if not warns else "WARNING"
    warn_code = "|".join(warns) if warns else "OK"

    with open(args.out, "w") as o:
        o.write(f"avg_depth\t{avg_depth:.2f}\n")
        o.write(f"total_pos\t{total}\n")
        o.write(f"low_dp\t{args.low_dp}\n")
        o.write(f"low_pos\t{low}\n")
        o.write(f"low_frac\t{low_frac:.6f}\n")
        o.write(f"zero_pos\t{zero}\n")
        o.write(f"zero_frac\t{zero_frac:.6f}\n")
        o.write(f"th_low_frac\t{args.low_frac:.6f}\n")
        o.write(f"th_zero_frac\t{args.zero_frac:.6f}\n")
        o.write(f"status\t{status}\n")
        o.write(f"warning\t{warn_code}\n")

if __name__ == "__main__":
    main()

