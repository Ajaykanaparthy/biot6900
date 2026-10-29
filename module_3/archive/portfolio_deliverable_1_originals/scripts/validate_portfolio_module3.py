"""Re-execute and independently audit the existing completed LUAD assignment.

Preserves original coursework before any update. Does not alter notebook code,
the Module 2 target file, Git state, or other modules. Run with course Python.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

import nbformat
from nbclient import NotebookClient
import numpy as np
import pandas as pd
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "validation"
ARCHIVE = ROOT / "archive" / "portfolio_deliverable_1_originals"
NOTEBOOK = ROOT / "BIOT6900_Module3_Assignment.ipynb"
CANONICAL = ROOT / "targets_luad_w3.csv"
TARGETS = ROOT.parent / "module_2" / "targets_luad.csv"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preserve(path):
    if not path.exists():
        return None
    relative = path.relative_to(ROOT)
    backup = ARCHIVE / relative
    if backup.exists() and backup.read_bytes() != path.read_bytes():
        backup = backup.with_name(f"{backup.stem}.{sha(path)[:12]}{backup.suffix}")
    backup.parent.mkdir(parents=True, exist_ok=True)
    if not backup.exists():
        shutil.copy2(path, backup)
    assert backup.read_bytes() == path.read_bytes()
    return backup.relative_to(ROOT).as_posix()


def save_table(name, table):
    path = OUT / name
    preserve(path)
    table.to_csv(path, index=False, lineterminator="\n")


def model():
    return make_pipeline(StandardScaler(), LogisticRegression(C=0.05, max_iter=5000))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    original_notebook_bytes = NOTEBOOK.read_bytes()
    original_csv_bytes = CANONICAL.read_bytes()
    original_csv = pd.read_csv(CANONICAL)
    original_sources = [cell.get("source", "") for cell in nbformat.read(NOTEBOOK, as_version=4).cells]
    backups = {p.name: preserve(p) for p in [NOTEBOOK, CANONICAL]}
    source_paths = [TARGETS, *sorted(DATA.glob("*"))]
    source_hashes = {str(p.relative_to(ROOT.parent)): sha(p) for p in source_paths if p.is_file()}

    # A fresh kernel and clear outputs implement Restart & Run All. All source
    # cells are retained verbatim, including the professor's TODO comments.
    notebook = nbformat.read(NOTEBOOK, as_version=4)
    for cell in notebook.cells:
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    print("Executing all cells in the existing completed notebook in a fresh kernel...", flush=True)
    client = NotebookClient(notebook, timeout=600, kernel_name="python3", allow_errors=False,
                            resources={"metadata": {"path": str(ROOT)}})
    try:
        client.execute(env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "MPLBACKEND": "Agg"})
    except Exception:
        # Keep completed canonical results even if a future environment fails.
        CANONICAL.write_bytes(original_csv_bytes)
        raise
    executed_sources = [cell.get("source", "") for cell in notebook.cells]
    assert original_sources == executed_sources, "Notebook source changed during execution"
    assert not any(output.get("output_type") == "error" for cell in notebook.cells
                   if cell.cell_type == "code" for output in cell.outputs)
    generated_csv_bytes = CANONICAL.read_bytes()
    reproduced = pd.read_csv(CANONICAL)
    # Preserve a separate exact regenerated artifact, then retain the historical
    # canonical CSV bytes. Do not silently replace a previously submitted result.
    preserve(OUT / "reproduced_targets_luad_w3.csv")
    (OUT / "reproduced_targets_luad_w3.csv").write_bytes(generated_csv_bytes)
    CANONICAL.write_bytes(original_csv_bytes)
    assert CANONICAL.read_bytes() == original_csv_bytes

    same_ranking = reproduced[["rank", "gene"]].equals(original_csv[["rank", "gene"]])
    same_neighbors = reproduced.nearest_drug_targets.equals(original_csv.nearest_drug_targets)
    deltas = {c: float(np.max(np.abs(reproduced[c] - original_csv[c])))
              for c in ["week2_score", "p_sm", "p_ab", "final"]}
    if not same_ranking or not same_neighbors or any(
        deltas[c] > (1e-4 if c in ("p_sm", "p_ab") else 1e-12) for c in deltas
    ):
        preserve(OUT / "reproduction_discrepancy.json")
        (OUT / "reproduction_discrepancy.json").write_text(json.dumps({
            "same_ranking": same_ranking, "same_neighbors": same_neighbors, "max_absolute_deltas": deltas,
            "canonical_original_preserved": True, "notebook_original_preserved": True,
        }, indent=2) + "\n")
        raise RuntimeError("Fresh notebook results differ materially; original notebook/CSV retained")

    # Independent calculations from the retained data pack, without importing
    # or executing any professor helper from the notebook.
    panel = pd.read_csv(DATA / "panel.tsv", sep="\t")
    proteome = pd.read_csv(DATA / "proteome.tsv", sep="\t")
    manifest = json.loads((DATA / "manifest.json").read_text())
    assert manifest["dry_run"] is False
    xp = np.load(DATA / "panel_emb_main.npy").astype("float32")
    xproteome = np.load(DATA / "proteome_emb_main.npy").astype("float32")
    assert xp.shape == (len(panel), 960) and xproteome.shape == (len(proteome), 960)
    assert np.isfinite(xp).all() and np.isfinite(xproteome).all()
    assert not proteome.gene.duplicated().any() and not panel.gene.duplicated().any()
    mine = pd.read_csv(TARGETS).sort_values("score", ascending=False).head(15).reset_index(drop=True)
    assert len(mine) == 15 and not mine.gene.duplicated().any()
    assert set(mine.gene) == set(original_csv.gene)
    lookup = {gene: i for i, gene in enumerate(proteome.gene)}
    assert all(gene in lookup for gene in mine.gene)
    rows = np.array([lookup[gene] for gene in mine.gene])
    mapped = proteome.iloc[rows].reset_index(drop=True)
    assert mapped.gene.tolist() == mine.gene.tolist()
    assert np.array_equal(mapped.length, mapped.sequence.str.len())
    keep = ~panel.gene.isin(mine.gene).to_numpy()
    train, xtrain = panel.loc[keep].reset_index(drop=True), xp[keep]
    overlaps = panel.loc[~keep, ["gene", "accession", "any_clinical", "family_group"]]
    assert len(overlaps) == 2 and set(overlaps.gene) == {"SCN7A", "GDF10"}
    assert not set(train.gene) & set(mine.gene)
    xtargets = xproteome[rows]
    y = train.any_clinical.to_numpy(dtype=int)
    groups = train.family_group.to_numpy()
    family_folds = GroupKFold(5)
    coverage = np.zeros(len(train), dtype=int)
    fold_rows = []
    for fold, (tr, te) in enumerate(family_folds.split(xtrain, y, groups), 1):
        assert not set(groups[tr]) & set(groups[te])
        coverage[te] += 1
        fold_rows.append({"fold": fold, "training_proteins": len(tr), "test_proteins": len(te),
                          "test_positive_labels": int(y[te].sum()), "overlapping_family_groups": 0})
    assert (coverage == 1).all()
    print("Independently reproducing random/family out-of-fold validation...", flush=True)
    random_p = cross_val_predict(model(), xtrain, y,
        cv=StratifiedKFold(5, shuffle=True, random_state=6900), method="predict_proba")[:, 1]
    family_p = cross_val_predict(model(), xtrain, y, cv=family_folds,
        groups=groups, method="predict_proba")[:, 1]
    metrics = pd.DataFrame([
        {"split": name, "AUROC": roc_auc_score(y, probabilities),
         "average_precision": average_precision_score(y, probabilities)}
        for name, probabilities in [("random", random_p), ("family", family_p)]
    ])
    family_metric = metrics.loc[metrics.split.eq("family")].iloc[0]
    assert round(family_metric.AUROC, 3) == 0.726
    assert round(family_metric.average_precision, 3) == 0.537
    psm = model().fit(xtrain, train.sm_clinical.astype(int)).predict_proba(xtargets)[:, 1]
    pab = model().fit(xtrain, train.ab_clinical.astype(int)).predict_proba(xtargets)[:, 1]
    clinical = train.any_clinical.to_numpy(dtype=bool)
    drugged = train.loc[clinical].reset_index(drop=True)
    xd = xtrain[clinical]
    similarities = ((xtargets / np.linalg.norm(xtargets, axis=1, keepdims=True)) @
                    (xd / np.linalg.norm(xd, axis=1, keepdims=True)).T)
    nearest_strings, neighbors = [], []
    for i, gene in enumerate(mine.gene):
        indices = np.argsort(-similarities[i])[:3]
        nearest_strings.append("; ".join(f"{drugged.iloc[j].gene} ({similarities[i, j]:.2f})" for j in indices))
        for rank, j in enumerate(indices, 1):
            neighbors.append({"gene": gene, "neighbor_rank": rank, "neighbor_gene": drugged.iloc[j].gene,
                              "neighbor_accession": drugged.iloc[j].accession,
                              "cosine_similarity": float(similarities[i, j])})
    result = mine[["gene", "score"]].rename(columns={"score": "module2_score"})
    result["module2_rank"] = np.arange(1, 16)
    result["accession"] = mapped.accession
    result["proteome_row"] = rows
    result["protein_length"] = mapped.length
    result["embedding_truncated"] = mapped.truncated_main
    result["p_sm"], result["p_ab"] = psm, pab
    result["nearest_drug_targets"] = nearest_strings
    result["week2_pct"] = result.module2_score.rank(pct=True)
    result["tractability_prior"] = result[["p_sm", "p_ab"]].max(axis=1)
    result["tractability_percentile"] = result.tractability_prior.rank(pct=True)
    result["final"] = 0.6 * result.week2_pct + 0.4 * result.tractability_percentile
    result["module3_rank"] = result.final.rank(ascending=False, method="first").astype(int)
    canonical = original_csv.set_index("gene")
    for _, row in result.iterrows():
        previous = canonical.loc[row.gene]
        assert previous["rank"] == row.module3_rank
        assert abs(previous.p_sm - row.p_sm) <= 1e-4 and abs(previous.p_ab - row.p_ab) <= 1e-4
        assert abs(previous.final - row.final) <= 1e-12
        assert previous.nearest_drug_targets == row.nearest_drug_targets

    anchor_table = pd.read_csv(DATA / "anchor_variants.tsv", sep="\t")
    anchor_rows = []
    for gene, variant in [("TP53", "R175H"), ("APOE", "C130R"), ("APOE", "R176C"), ("TREM2", "R47H")]:
        source = pd.read_csv(DATA / f"variant_scores_{gene}.tsv", sep="\t")
        one = source.loc[source.variant.eq(variant)]
        assert len(one) == 1
        actual = one.iloc[0]
        declared = anchor_table.loc[anchor_table.gene.eq(gene) & anchor_table.variant.eq(variant)].iloc[0]
        percentile = float(100 * source.score.le(actual.score).mean())
        assert actual.score == declared.score and round(percentile, 1) == declared.percentile_in_protein
        precursor = proteome.loc[proteome.gene.eq(gene)].iloc[0]
        assert precursor.sequence[int(actual.position) - 1] == actual.wt
        anchor_rows.append({"gene": gene, "variant": variant, "score": float(actual.score),
                            "percentile_source": float(declared.percentile_in_protein),
                            "percentile_recomputed": percentile, "matching_source_score": True,
                            "wild_type_residue_matches_precursor_sequence": True,
                            "source_substitution_count": len(source)})

    save_table("target_mapping_and_scores.csv", result.sort_values("module3_rank"))
    save_table("nearest_drug_targets.csv", pd.DataFrame(neighbors))
    save_table("anchor_variant_audit.csv", pd.DataFrame(anchor_rows))
    save_table("validation_metrics.csv", metrics)
    save_table("family_fold_audit.csv", pd.DataFrame(fold_rows))
    save_table("training_overlap_removed.csv", overlaps)
    save_table("out_of_fold_predictions.csv", pd.DataFrame({"gene": train.gene, "family_group": groups,
        "any_clinical": y, "random_oof_probability": random_p, "family_oof_probability": family_p}))
    assert all(sha(ROOT.parent / p) == value for p, value in source_hashes.items()), "Source inputs changed"
    assert CANONICAL.read_bytes() == original_csv_bytes
    # Save fresh outputs only after successful comparison and independent checks.
    preserve(NOTEBOOK)
    nbformat.write(notebook, NOTEBOOK)
    assert [cell.get("source", "") for cell in nbformat.read(NOTEBOOK, as_version=4).cells] == original_sources
    report = {
        "status": "PASS", "notebook": NOTEBOOK.name, "fresh_kernel_restart_run_all": True,
        "executed_code_cells": sum(c.cell_type == "code" for c in notebook.cells),
        "notebook_source_cells_unchanged": True,
        "original_notebook_sha256": hashlib.sha256(original_notebook_bytes).hexdigest(),
        "executed_notebook_sha256": sha(NOTEBOOK), "backups": backups,
        "canonical_original_csv_bytes_preserved": True, "canonical_csv_sha256": sha(CANONICAL),
        "fresh_csv_byte_identical": generated_csv_bytes == original_csv_bytes,
        "same_15_genes_and_ranks": same_ranking, "same_nearest_neighbor_strings": same_neighbors,
        "maximum_numeric_differences_from_original_csv": deltas,
        "source_data_hashes": source_hashes,
        "mapped_targets": 15, "unmapped_targets": [], "training_panel_before": len(panel),
        "training_panel_after": len(train), "targets_removed_from_training": overlaps.gene.tolist(),
        "five_family_folds_have_no_group_overlap": True,
        "validation_metrics": metrics.to_dict(orient="records"),
        "chance_average_precision": float(y.mean()),
        "ranking_method": "0.6*percentile_rank(Module2_score) + 0.4*percentile_rank(max(p_sm,p_ab)); ascending percentiles; final descending rank(method='first')",
        "ranking_tie_handling": "Retain professor's existing arithmetic and method='first'; tiny floating point differences can distinguish apparent printed ties",
        "embedding_model": manifest["embedding_models"]["main"], "dry_run": manifest["dry_run"],
        "truncated_target_embeddings": mapped.loc[mapped.truncated_main, "gene"].tolist(),
        "anchor_variants_verified": anchor_rows,
        "limitations": ["Course-supplied embeddings and labels are reused; model inference was not rerun from sequences",
                        "Clinical-precedence priors are not calibrated probabilities of LUAD drug success",
                        "Group validation applies to training labels, not independent LUAD therapeutic efficacy",
                        "COL6A6 and TNXB embeddings are flagged truncated in the supplied data",
                        "Variant log-likelihood scores and percentiles do not establish pathogenicity"],
        "versions": {"python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__,
                     "scikit_learn": sklearn.__version__, "nbformat": nbformat.__version__},
    }
    preserve(OUT / "validation_report.json")
    (OUT / "validation_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in ["status", "fresh_kernel_restart_run_all", "executed_code_cells",
        "canonical_original_csv_bytes_preserved", "fresh_csv_byte_identical", "maximum_numeric_differences_from_original_csv",
        "targets_removed_from_training", "validation_metrics", "truncated_target_embeddings"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
