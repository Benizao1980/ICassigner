import subprocess
import sys
from pathlib import Path

def test_icassigner_smoke(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]

    tree = repo_root / "example" / "Acinetobacter-coreML.nwk"
    meta = repo_root / "example" / "metadata.csv"
    script = repo_root / "ICassigner.py"

    outdir = tmp_path / "outputs"

    cmd = [
        sys.executable,
        str(script),
        "--tree", str(tree),
        "--metadata", str(meta),
        "--tip_col", "sample_id",
        "--ic_col", "IC",
        "--outdir", str(outdir)
    ]

    result = subprocess.run(cmd, capture_output=True)

    assert result.returncode == 0, result.stderr.decode()
    assert (outdir / "ICassigner_summary.txt").exists()
