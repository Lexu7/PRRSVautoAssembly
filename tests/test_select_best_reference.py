import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "select_best_reference.py"


def test_selects_largest_cumulative_bitscore_and_breaks_ties_by_id(tmp_path):
    blast = tmp_path / "hits.tsv"
    blast.write_text(
        "contig1\tref_b\t99\t100\t100\t15000\t1\t100\t1\t100\t1e-20\t50\n"
        "contig2\tref_b\t99\t100\t100\t15000\t1\t100\t1\t100\t1e-20\t25\n"
        "contig1\tref_a\t99\t100\t100\t15000\t1\t100\t1\t100\t1e-20\t75\n"
    )
    selected = tmp_path / "selected.txt"
    scores = tmp_path / "scores.tsv"

    subprocess.run(
        [sys.executable, str(SCRIPT), "--blast", str(blast), "--out", str(selected), "--scores", str(scores)],
        check=True,
    )

    assert selected.read_text() == "ref_a\n"
    assert scores.read_text().splitlines()[1] == "ref_a\t75.000000"
