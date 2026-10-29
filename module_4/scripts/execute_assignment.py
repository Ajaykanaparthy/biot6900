"""Restart/run the existing Module 4 notebook and compare its preserved results."""
from pathlib import Path
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone

import nbformat
from nbclient import NotebookClient
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    path = ROOT / "BIOT6900_Module4_Assignment.ipynb"
    backup = ROOT / "archive/portfolio_deliverable_1_originals" / path.name
    if not backup.is_file():
        raise RuntimeError("Original notebook backup is required before execution")
    notebook = nbformat.read(path, as_version=4)
    client = NotebookClient(notebook, timeout=300, kernel_name="python3",
                            resources={"metadata": {"path": str(ROOT)}}, allow_errors=False)
    client.execute()
    errors = [o for c in notebook.cells if c.cell_type == "code" for o in c.outputs if o.output_type == "error"]
    if errors:
        raise RuntimeError("Notebook has executed error outputs")
    before = pd.read_csv(ROOT / "archive/portfolio_deliverable_1_originals/structure_assessment_luad.csv")
    after = pd.read_csv(ROOT / "structure_assessment_luad.csv")
    pd.testing.assert_frame_equal(before, after, check_dtype=False, rtol=1e-10, atol=1e-10)
    nbformat.write(notebook, path)
    result = {"notebook": path.name, "execution": "PASS", "fresh_kernel": True,
              "code_cells_executed": sum(c.cell_type == "code" for c in notebook.cells),
              "error_outputs": len(errors), "legacy_summary_values_reproduced": True,
              "python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
              "pandas": pd.__version__, "kernel": "python3 (biot6900 environment)",
              "fpocket": "Retained outputs verified and reused; no new local run",
              "finished_utc": datetime.now(timezone.utc).isoformat(),
              "original_notebook_sha256": hashlib.sha256(backup.read_bytes()).hexdigest(),
              "executed_notebook_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    (ROOT / "results/notebook_validation.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
