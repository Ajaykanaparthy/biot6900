"""Join verified existing Modules 2-4 outputs without changing their rankings.

Run after each module's independent validation has finished. Missing structural
assessments remain blank and explicitly marked unassessed. No values are imputed.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "portfolio_deliverable_1"
INPUTS = {
    "module2_targets": ROOT / "module_2/targets_luad.csv",
    "module2_integrated": ROOT / "module_2/data/luad_integrated_gene_level.csv",
    "module3_ranking": ROOT / "module_3/targets_luad_w3.csv",
    "module3_validation": ROOT / "module_3/validation/target_mapping_and_scores.csv",
    "module4_structure": ROOT / "module_4/results/verified_structure_summary.csv",
}


def flag(value) -> bool:
    if pd.isna(value):
        raise ValueError("Required evidence flag is missing; it must not be imputed")
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if str(value).lower() in {"true", "1"}:
        return True
    if str(value).lower() in {"false", "0"}:
        return False
    raise ValueError(f"Unrecognized boolean: {value!r}")


def new_text(path: Path, text: str) -> None:
    if path.exists():
        if path.read_text(encoding="utf-8") == text:
            return
        raise FileExistsError(f"Refusing to overwrite an existing different output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        handle.write(text)


def main() -> None:
    data = {name: pd.read_csv(path) for name, path in INPUTS.items()}
    targets = data["module2_targets"]
    integrated = data["module2_integrated"].set_index("gene", verify_integrity=True)
    ranking = data["module3_ranking"]
    validation = data["module3_validation"].set_index("gene", verify_integrity=True)
    structures = data["module4_structure"].set_index("gene", verify_integrity=True)
    assert len(targets) == len(ranking) == len(validation) == 15
    assert len(integrated) == 6758
    assert set(targets.gene) == set(ranking.gene) == set(validation.index)
    assert list(ranking["rank"]) == list(range(1, 16))
    assert set(structures.index) == {"SCN7A", "MASP1"}
    assert int(integrated["concordant"].map(flag).sum()) == 4918
    assert len(integrated) - int(integrated["concordant"].map(flag).sum()) == 1840
    target_ranks = {gene: i + 1 for i, gene in enumerate(targets.gene)}
    rows = []
    for rank_row in ranking.to_dict("records"):
        gene = rank_row["gene"]
        m2 = integrated.loc[gene]
        m3 = validation.loc[gene]
        assert int(m2.module2_rank) == target_ranks[gene]
        assert int(m3.module3_rank) == int(rank_row["rank"])
        assert np.isclose(float(m2.score), float(rank_row["week2_score"]), rtol=1e-12, atol=1e-12)
        for score in ("p_sm", "p_ab", "final"):
            assert np.isclose(float(m3[score]), float(rank_row[score]), rtol=2e-5, atol=2e-7)
        concordant = flag(m2.concordant)
        truncated = flag(m3.embedding_truncated)
        row = {
            "gene": gene,
            "module2_rank": target_ranks[gene],
            "rna_log2fc": m2.rna_lfc,
            "rna_p": m2.rna_p,
            "rna_q_bh": m2.rna_q_bh,
            "rna_valid_pairs": m2.rna_n_pairs,
            "rna_positive_evidence": flag(m2.evidence_rna),
            "protein_log2fc": m2.prot_lfc,
            "protein_p": m2.prot_p,
            "protein_q_bh": m2.prot_q_bh,
            "protein_valid_pairs": m2.prot_n_pairs,
            "protein_positive_evidence": flag(m2.evidence_protein),
            "mutation_frequency": m2.mutation_freq,
            "genomic_positive_evidence": flag(m2.evidence_genomic),
            "rna_protein_concordant": concordant,
            "module2_score": m2.score,
            "module3_rank": int(rank_row["rank"]),
            "uniprot": m3.accession,
            "protein_length": m3.protein_length,
            "embedding_truncated": truncated,
            "small_molecule_tractability_score": rank_row["p_sm"],
            "antibody_tractability_score": rank_row["p_ab"],
            "nearest_known_drug_target_similarities": rank_row["nearest_drug_targets"],
            "module3_final_score": rank_row["final"],
            "structure_assessed": gene in structures.index,
            "mean_pLDDT": None,
            "percent_residues_below_pLDDT50": None,
            "experimental_structure": None,
            "confident_region_rmsd_angstrom": None,
            "best_pocket_id": None,
            "fpocket_score": None,
            "fpocket_druggability_score": None,
            "pocket_mean_pLDDT": None,
            "pocket_within_PAE_angstrom": None,
            "pocket_confidence": "not assessed",
            "structure_druggability_call": "not assessed",
        }
        limitations = [
            "Gene-level association and tractability predictions do not establish LUAD dependency or therapeutic benefit",
            "Mutation recurrence is not a driver-significance test",
        ]
        if truncated:
            limitations.append("Course main protein embedding is truncated; it does not represent the complete sequence")
        if not concordant:
            limitations.append("RNA and protein directions disagree; matched orthogonal validation is needed")
        if int(m2.prot_n_pairs) < 10:
            limitations.append(f"Protein effect rests on only {int(m2.prot_n_pairs)} valid pairs and fails the portfolio evidence support threshold")
        if gene in structures.index:
            structure = structures.loc[gene]
            assert str(structure.uniprot) == str(m3.accession)
            required_structure_numbers = [
                structure.mean_pLDDT, structure.pct_below_50, structure.rmsd_confident,
                structure.best_pocket_id, structure.best_fpocket_score,
                structure.best_pocket_druggability_score,
                structure.best_pocket_mean_pLDDT, structure.best_pocket_PAE_within,
            ]
            assert np.isfinite(np.asarray(required_structure_numbers, dtype=float)).all()
            row.update({
                "mean_pLDDT": structure.mean_pLDDT,
                "percent_residues_below_pLDDT50": structure.pct_below_50,
                "experimental_structure": structure.pdb_id,
                "confident_region_rmsd_angstrom": structure.rmsd_confident,
                "best_pocket_id": int(structure.best_pocket_id),
                "fpocket_score": structure.best_fpocket_score,
                "fpocket_druggability_score": structure.best_pocket_druggability_score,
                "pocket_mean_pLDDT": structure.best_pocket_mean_pLDDT,
                "pocket_within_PAE_angstrom": structure.best_pocket_PAE_within,
                "pocket_confidence": structure.best_pocket_verdict,
                "structure_druggability_call": structure.druggability_call,
            })
            if not pd.isna(structure.limitations):
                limitations.append(str(structure.limitations))
        else:
            limitations.append("No Module 4 structural or pocket assessment was performed for this candidate")
        if gene == "COL6A6":
            rationale = "Preserve Module 3 rank 1 and Module 2 rank 3; concordant decreases support prioritization, while modest tractability scores, truncation, and absent structural validation limit a drug-target claim"
        elif gene == "SCN7A":
            rationale = "Preserve Module 3 rank 3 and Module 2 rank 5; relatively strong small-molecule classifier score and a locally confident pocket support follow-up, subject to uncertain inter-domain placement"
        elif gene == "MASP1":
            rationale = "Preserve Module 3 rank 8 and Module 2 rank 11; both modality scores and a retained pocket support follow-up, while the experimental comparison covers a mutant catalytic region"
        elif gene == "ZNF536":
            rationale = "Retain the notable discordant Module 2 rank-2 hit; low modality scores and protein evidence limitations argue for expression validation before therapeutic prioritization"
        else:
            rationale = "Preserve the documented Module 3 rank combining Module 2 evidence and modality tractability; structural and disease-functional validation remain outstanding"
        row["final_prioritization_rationale"] = rationale
        row["limitations"] = "; ".join(limitations)
        rows.append(row)
    output = pd.DataFrame(rows)
    assert output.gene.nunique() == 15
    for gene, expected in {"COL6A6": (3, 1), "SCN7A": (5, 3), "MASP1": (11, 8)}.items():
        row = output.set_index("gene").loc[gene]
        assert (int(row.module2_rank), int(row.module3_rank)) == expected
    assert int(output.structure_assessed.sum()) == 2
    assert output.loc[~output.structure_assessed, "mean_pLDDT"].isna().all()
    assert output.loc[~output.structure_assessed, "fpocket_druggability_score"].isna().all()
    csv_text = output.to_csv(index=False, lineterminator="\n", float_format="%.12g")
    lines = [
        "# Verified cross-module LUAD target summary",
        "",
        "Rows follow the existing Module 3 ranking. No new composite ranking was invented.",
        "RNA and protein values are paired tumor-minus-normal log2 effects; mutation frequency is a cohort proportion.",
        "Tractability outputs are course classifier scores, not validated probabilities of clinical success.",
        "Displayed classifier scores are rounded to three decimals; exact nonzero values are retained in the CSV.",
        "A dash means the structure/pocket was not assessed, rather than a score of zero.",
        "",
        "| Gene | M2 rank | M3 rank | RNA log2FC | Protein log2FC | Mutation % | Concordant | SM score | Antibody score | Mean pLDDT | Pocket druggability / confidence |",
        "|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---|",
    ]
    for row in rows:
        mean = "-" if row["mean_pLDDT"] is None else f"{row['mean_pLDDT']:.2f}"
        pocket = "-" if row["fpocket_druggability_score"] is None else f"{row['fpocket_druggability_score']:.3f} / {row['pocket_confidence']}"
        lines.append(
            f"| {row['gene']} | {row['module2_rank']} | {row['module3_rank']} | {row['rna_log2fc']:.3f} | "
            f"{row['protein_log2fc']:.3f} | {100 * row['mutation_frequency']:.2f} | "
            f"{'yes' if row['rna_protein_concordant'] else 'no'} | {row['small_molecule_tractability_score']:.3f} | "
            f"{row['antibody_tractability_score']:.3f} | {mean} | {pocket} |"
        )
    lines += [
        "",
        "The [CSV](target_prioritization_summary.csv) includes all pair counts, p/q-values, thresholded evidence flags, UniProt mappings, embedding truncation flags, separate fpocket Score and Druggability Score, and candidate-specific rationales/limitations.",
        "",
        "A high Module 2 rank measures relative dysregulation and mutation recurrence. It does not demonstrate ligand binding, tumor dependency, benefit from inhibiting a downregulated gene, or an acceptable safety margin.",
        "",
    ]
    provenance = {
        "purpose": "Synthesis of verified existing Modules 2-4; no new target ranking",
        "targets": len(output),
        "structure_assessed": int(output.structure_assessed.sum()),
        "missing_structure_values_imputed": False,
        "input_sha256": {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in INPUTS.values()},
        "csv_sha256": hashlib.sha256(csv_text.encode("utf-8")).hexdigest(),
    }
    new_text(DEST / "target_prioritization_summary.csv", csv_text)
    new_text(DEST / "target_prioritization_summary.md", "\n".join(lines))
    new_text(DEST / "summary_provenance.json", json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"targets": len(output), "structure_assessed": 2,
                      "module2_joined": len(integrated), "concordant": 4918,
                      "module3_order_preserved": True, "output": "portfolio_deliverable_1/target_prioritization_summary.csv"}, indent=2))


if __name__ == "__main__":
    main()
