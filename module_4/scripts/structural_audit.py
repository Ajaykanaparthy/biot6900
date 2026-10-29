"""Verify retained structural inputs and create reproducible Module 4 figures.

Existing fpocket output is parsed, never represented as a fresh fpocket run.
The original command/version and Explorer job log are unavailable.
"""
from pathlib import Path
import hashlib
import json
import re
import zipfile

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ensure_local_inputs():
    """Restore only missing, preserved inputs from the accompanying ZIP."""
    archive = ROOT / "reproducibility/luad_structure_inputs.zip"
    if not archive.exists():
        return
    with zipfile.ZipFile(archive) as handle:
        for entry in handle.infolist():
            if entry.is_dir():
                continue
            target = (ROOT / entry.filename).resolve()
            if not target.is_relative_to((ROOT / "structures").resolve()):
                raise ValueError("Archive entry is outside the structures directory")
            payload = handle.read(entry)
            if target.exists():
                if target.read_bytes() != payload:
                    raise ValueError(f"Retained input differs from archived input: {target}")
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open("xb") as stream:
                    stream.write(payload)


def _atoms(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM") and line[16] in (" ", "A"):
            key = (line[21], int(line[22:26]), line[26].strip(), line[12:16].strip())
            if key in result:
                raise ValueError(f"Duplicate atom identity: {path} {key}")
            result[key] = (line[17:20].strip(), np.array([
                float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    return result


def validate_fpocket_cache(pdb_path, out_dir):
    """Check every retained lining atom against its exact source AF model."""
    pdb_path, out_dir = Path(pdb_path).resolve(), Path(out_dir).resolve()
    info_path = out_dir / f"{pdb_path.stem}_info.txt"
    text = info_path.read_text()
    pocket_ids = [int(x) for x in re.findall(r"(?m)^Pocket\s+(\d+)\s*:", text)]
    if not pocket_ids or len(set(pocket_ids)) != len(pocket_ids):
        raise ValueError(f"Missing/duplicate pocket IDs in {info_path}")
    source = _atoms(pdb_path)
    checked, maximum = 0, 0.0
    hashes = {pdb_path.relative_to(ROOT).as_posix(): sha256(pdb_path),
              info_path.relative_to(ROOT).as_posix(): sha256(info_path)}
    for pocket_id in pocket_ids:
        atom_path = out_dir / "pockets" / f"pocket{pocket_id}_atm.pdb"
        vertex_path = out_dir / "pockets" / f"pocket{pocket_id}_vert.pqr"
        atoms = _atoms(atom_path)
        if not atoms or not vertex_path.exists():
            raise ValueError(f"Incomplete fpocket output for pocket {pocket_id}")
        for key, (resname, xyz) in atoms.items():
            if key not in source or source[key][0] != resname:
                raise ValueError(f"Pocket atom identity differs from AF source: {atom_path} {key}")
            delta = float(np.linalg.norm(xyz - source[key][1]))
            if delta > 0.002:
                raise ValueError(f"Pocket atom coordinates differ from source: {atom_path} {key}")
            maximum = max(maximum, delta)
            checked += 1
        for path in (atom_path, vertex_path):
            hashes[path.relative_to(ROOT).as_posix()] = sha256(path)
    return {"source": "retained fpocket-format outputs attributed by the student to Explorer",
            "fresh_fpocket_execution": False,
            "original_job_log_command_and_version_available": False,
            "pockets": len(pocket_ids), "lining_atom_entries_checked": checked,
            "maximum_coordinate_difference_A": maximum,
            "all_lining_atom_identities_and_coordinates_match_AF_source": True,
            "sha256": hashes}


def pocket_quality(lining_residues, ca_table, pae):
    """Use residue-weighted pLDDT and all ordered PAE pairs, including diagonal."""
    by_residue = ca_table.set_index("resnum").b
    ids = np.asarray(lining_residues, dtype=int)
    if not len(ids) or len(set(ids)) != len(ids) or not set(ids).issubset(by_residue.index):
        raise ValueError("Pocket lining residues are empty, duplicated, or outside the model")
    # The retained AF models have contiguous UniProt numbering, checked explicitly.
    if list(ca_table.resnum) != list(range(1, len(ca_table) + 1)):
        raise ValueError("PAE lookup requires contiguous source residue numbering")
    mean = float(by_residue.loc[ids].mean())
    within = float(pae[np.ix_(ids - 1, ids - 1)].mean())
    verdict = "reject" if mean < 70 else "suspect" if within > 10 else "keep"
    return pd.Series({"mean_pLDDT": mean, "PAE_within": within, "verdict": verdict})


def save_structural_figures(targets, ca, pae, all_pockets):
    directory = ROOT / "Figures"
    directory.mkdir(exist_ok=True)
    cmap = ListedColormap(["#FF7D45", "#FFDB13", "#65CBF3", "#0053D6"])
    norm = BoundaryNorm([0, 50, 70, 90, 100.01], cmap.N)
    captions = ["# Structural figure captions", ""]
    for gene, accession in targets.items():
        table = ca[gene]
        pockets = all_pockets[gene].sort_values("Druggability Score", ascending=False)
        kept = pockets[pockets.verdict == "keep"]
        if kept.empty:
            raise ValueError(f"No kept pocket available for the structural figure: {gene}")
        best = kept.iloc[0]
        xyz = table[["x", "y", "z"]].to_numpy()
        residue_mask = table.resnum.isin(best.lining_residues).to_numpy()
        fig = plt.figure(figsize=(14, 10), constrained_layout=True)
        grid = fig.add_gridspec(2, 2)
        ax3d = fig.add_subplot(grid[0, 0], projection="3d")
        segments = np.stack([xyz[:-1], xyz[1:]], axis=1)
        lines = Line3DCollection(segments, cmap=cmap, norm=norm, linewidth=1.1)
        lines.set_array(table.b.to_numpy()[:-1])
        ax3d.add_collection3d(lines)
        ax3d.scatter(*xyz[residue_mask].T, color="#BB1C83", s=32, depthshade=False)
        lower, upper = xyz.min(axis=0), xyz.max(axis=0)
        for setter, lo, hi in zip((ax3d.set_xlim, ax3d.set_ylim, ax3d.set_zlim), lower, upper):
            setter(lo - 3, hi + 3)
        ax3d.set_box_aspect(np.maximum(upper - lower, 1))
        ax3d.view_init(elev=18, azim=35)
        ax3d.set_axis_off()
        ax3d.set_title(f"A  AlphaFold C-alpha trace\nMagenta: lining residues of pocket {int(best.pocket)}", loc="left")

        ax_conf = fig.add_subplot(grid[0, 1])
        ax_conf.scatter(table.resnum, table.b, c=table.b, cmap=cmap, norm=norm, s=4, rasterized=True)
        ax_conf.scatter(table.resnum[residue_mask], table.b[residue_mask], color="#BB1C83", s=25,
                        label=f"Pocket {int(best.pocket)} lining residues", zorder=3)
        ax_conf.axhline(70, ls="--", lw=1, color="#333333", label="pLDDT 70")
        ax_conf.set(xlim=(0, len(table) + 1), ylim=(0, 103), xlabel="UniProt residue", ylabel="pLDDT")
        ax_conf.set_title("B  Local model confidence", loc="left")
        ax_conf.legend(frameon=False, loc="lower center", fontsize=9)

        ax_pae = fig.add_subplot(grid[1, 0])
        n = len(table)
        im = ax_pae.imshow(pae[gene], origin="upper", extent=(0.5, n + 0.5, n + 0.5, 0.5),
                           cmap="viridis", vmin=0, vmax=30, interpolation="nearest")
        ax_pae.set(xlabel="Aligned-on residue", ylabel="Predicted residue")
        ax_pae.set_title("C  Predicted aligned error", loc="left")
        fig.colorbar(im, ax=ax_pae, shrink=0.8, label="PAE (Angstrom)")

        ax_pk = fig.add_subplot(grid[1, 1])
        top = pockets.head(5).iloc[::-1]
        colors = ["#B7443D" if v == "reject" else "#D9A021" if v == "suspect" else "#20746B" for v in top.verdict]
        ax_pk.barh(range(len(top)), top["Druggability Score"], color=colors)
        ax_pk.set_yticks(range(len(top)), [f"Pocket {int(p)}" for p in top.pocket])
        ax_pk.set(xlim=(0, 1.42), xlabel="fpocket druggability score")
        ax_pk.set_title("D  Top five pockets and confidence verdicts", loc="left")
        for y, (_, row) in enumerate(top.iterrows()):
            ax_pk.text(float(row["Druggability Score"]) + 0.02, y,
                       f"{row['Druggability Score']:.3f}  {row.verdict}\npLDDT {row.mean_pLDDT:.1f}; PAE {row.PAE_within:.2f}",
                       va="center", fontsize=9)
        legend = [Line2D([0], [0], marker="o", linestyle="", color=cmap(i), label=label)
                  for i, label in enumerate(("pLDDT <50", "50-70", "70-90", ">90"))]
        ax3d.legend(handles=legend, loc="lower left", ncol=2, frameon=False, fontsize=9)
        fig.suptitle(f"{gene} ({accession}): structural confidence and retained fpocket evidence", fontsize=16)
        name = f"{gene.lower()}_structure_confidence_pockets.png"
        fig.savefig(directory / name, dpi=200, bbox_inches="tight")
        plt.close(fig)
        captions += [f"## {name}", "",
            f"**{gene}.** A: actual C-alpha coordinates colored by source AlphaFold pLDDT; "
            f"magenta marks the lining residues of the highest-druggability kept pocket ({int(best.pocket)}). "
            "This is a backbone trace, not a molecular-surface rendering. B: local confidence along "
            "the sequence, with the same pocket highlighted. C: full source PAE matrix. D: five "
            "highest-druggability pockets, retaining rejected candidates. Confidence follows the "
            "professor's rules: reject mean lining pLDDT <70; otherwise suspect mean lining PAE >10 Angstrom; "
            "otherwise keep. Mean lining PAE includes both directions and diagonal entries, as in the lab helper. "
            "A geometry/confidence-supported pocket does not demonstrate ligand binding, selectivity, or LUAD efficacy.", ""]
    (directory / "captions.md").write_text("\n".join(captions), encoding="utf-8")


def export_verified_tables(targets, target_info, ca, plddt_summary, domains, inter_pae,
                           comparison, all_pockets, calls, cache_audits, ligand_pockets):
    results = ROOT / "results"
    results.mkdir(exist_ok=True)
    rows, all_rows = [], []
    for gene, accession in targets.items():
        pockets = all_pockets[gene].sort_values("Druggability Score", ascending=False)
        kept = pockets[pockets.verdict == "keep"]
        if kept.empty:
            raise ValueError(f"No kept pocket for {gene}; update the explicit assessment")
        best = kept.iloc[0]
        for _, p in pockets.iterrows():
            row = p.to_dict()
            row["gene"] = gene
            row["lining_residues"] = ";".join(str(r) for r in row["lining_residues"])
            all_rows.append(row)
        confidence = plddt_summary.set_index("gene").loc[gene]
        limitation = (
            "High inter-component PAE limits relative placement; membrane environment and ligand binding not validated; only one confident component has sufficient experimental pairs"
            if gene == "SCN7A" else
            "Best pocket is close to pLDDT threshold; 4IGD is R448Q zymogen fragment 298-699, not full-length wild type; no retained functional inhibitor validates this pocket")
        rows.append({"gene": gene, "uniprot": accession,
            "module2_rank": target_info[gene]["module2_rank"],
            "mean_pLDDT": confidence.mean_pLDDT, "pct_below_50": confidence.pct_below_50,
            "residues": len(ca[gene]), "n_domains": len(domains[gene]),
            "inter_domain_PAE": inter_pae[gene], "pdb_id": comparison[gene]["pdb_id"],
            "pdb_chain": "A", "rmsd_confident": comparison[gene]["rmsd_confident"],
            "pairs_used": comparison[gene]["pairs_used"], "best_pocket_id": int(best.pocket),
            "best_fpocket_score": best["Score"],
            "best_pocket_druggability_score": best["Druggability Score"],
            "best_pocket_mean_pLDDT": best.mean_pLDDT, "best_pocket_PAE_within": best.PAE_within,
            "best_pocket_verdict": best.verdict, "druggability_call": calls[gene],
            "limitations": limitation})
        domains[gene].to_csv(results / f"{gene.lower()}_confidence_components.csv", index=False)
        if gene in ligand_pockets:
            ligand_pockets[gene].to_csv(results / f"{gene.lower()}_ligand_site_check.csv", index=False)
    pd.DataFrame(rows).to_csv(results / "verified_structure_summary.csv", index=False)
    pd.DataFrame(all_rows).to_csv(results / "pocket_quality_luad.csv", index=False)
    pd.DataFrame(all_rows).groupby("gene", sort=False).head(5).to_csv(results / "top_five_pockets_luad.csv", index=False)
    plddt_summary.to_csv(results / "plddt_summary_luad.csv", index=False)
    manifest = {"scope": "Audit and reuse of existing structural evidence; no fabricated fpocket run",
        "confidence_rule": "reject if mean_pLDDT <70; else suspect if PAE_within >10 Angstrom; else keep",
        "mean_PAE_definition": "All ordered pairs of unique lining residues, including diagonal",
        "cache_audits": cache_audits,
        "raw_input_sha256": {p.relative_to(ROOT).as_posix(): sha256(p) for p in (ROOT / "structures").rglob("*")
                             if p.is_file() and (p.name.startswith(("AF-Q01118", "AF-P48740")) or any(
                                 s in p.parts for s in ("AF-Q01118_out", "AF-P48740_out")) or p.name in ("7TJ8.pdb", "4IGD.pdb", "1TUP.pdb", "AF-P04637.pdb"))},
        "notes": ["Original Explorer job log, fpocket version, and command are not available",
                  "Exact atom identity/coordinate matches verify the cached pockets refer to the preserved AF inputs",
                  "find_domains is the provided connected-component heuristic, not an experimentally established domain annotation",
                  "4IGD is an engineered MASP1 catalytic-region zymogen fragment (298-699; R448Q)",
                  "The full published biological assembly and a drug-bound physiological state are not modeled"]}
    (results / "structural_input_provenance.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def audit_p53_anchor(read_pdb, ca_table, kabsch):
    """Reproduce the existing Part 2 p53 anchor from its retained structures."""
    af = ca_table(read_pdb(ROOT / "structures/AF-P04637.pdb")).set_index("resnum")
    experimental = read_pdb(ROOT / "structures/1TUP.pdb")
    chain = experimental[(experimental.chain == "B") & (experimental.record == "ATOM")]
    xca = ca_table(chain).set_index("resnum")
    common = sorted(set(af.index) & set(xca.index))
    confident = [r for r in common if af.loc[r, "b"] >= 70]
    _, _, rmsd = kabsch(af.loc[confident, ["x", "y", "z"]].to_numpy(),
                        xca.loc[confident, ["x", "y", "z"]].to_numpy())
    zinc = experimental[experimental.resname == "ZN"]
    dna = experimental[experimental.resname.isin(["DA", "DC", "DG", "DT", "A", "C", "G", "T"])]
    def distance(atoms, other):
        if atoms.empty or other.empty:
            return np.nan
        delta = atoms[["x", "y", "z"]].to_numpy()[:, None] - other[["x", "y", "z"]].to_numpy()[None]
        return float(np.sqrt((delta ** 2).sum(-1)).min())
    rows = []
    for residue in [175, 176, 179, 220, 238, 242, 245, 248, 249, 273, 282]:
        side = chain[(chain.resnum == residue) & ~chain.atom.isin(["N", "C", "O"])]
        d_zn, d_dna = distance(side, zinc), distance(side, dna)
        role = "zinc ligand" if d_zn < 3 else "DNA contact" if d_dna < 4.5 else "near DNA" if d_dna < 8 else "structural / other"
        rows.append({"residue": residue, "amino_acid": side.resname.iloc[0],
                     "pLDDT": float(af.loc[residue, "b"]), "distance_to_zinc_A": d_zn,
                     "distance_to_DNA_A": d_dna, "role": role})
    pd.DataFrame(rows).to_csv(ROOT / "results/p53_anchor_verified.csv", index=False)
    stats = {"pdb": "1TUP", "chain": "B", "common_residues": len(common),
             "confident_residues": len(confident), "rmsd_confident_A": float(rmsd)}
    (ROOT / "results/p53_anchor_summary.json").write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    return stats
