"""Recover the existing LUAD gene-level analysis and portfolio figures.

Uses the original five local CPTAC-LUAD matrices and original equal-rank scoring.
The first run retrieves genuine MyGene annotations; subsequent runs use the
retained annotation response. Existing target CSVs are validated, never replaced.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
import requests
from scipy import stats

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "luad_raw"
DATA = ROOT / "data"
FIGURES = ROOT / "Figures"
MAPPING = DATA / "luad_ensembl_symbol_mapping.tsv"
RESPONSE = DATA / "luad_mygene_response.json"
FALLBACK_RESPONSE = DATA / "luad_ensembl_fallback_response.json"

# Prespecified before examining set intersections. These do NOT alter ranking.
EVIDENCE_THRESHOLDS = {
    "rna": {"max_bh_q": 0.05, "min_abs_log2fc": 1.0, "min_valid_pairs": 10},
    "protein": {"max_bh_q": 0.05, "min_abs_log2fc": 1.0, "min_valid_pairs": 10},
    "genomic": {"min_mutation_frequency": 0.05},
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bh_qvalues(values):
    """Benjamini-Hochberg over finite tests only; missing p-values remain missing."""
    values = np.asarray(values, dtype=float)
    out = np.full(len(values), np.nan)
    positions = np.flatnonzero(np.isfinite(values))
    ordered = positions[np.argsort(values[positions], kind="stable")]
    n = len(ordered)
    adjusted = values[ordered] * n / np.arange(1, n + 1)
    out[ordered] = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
    return out


def raw_summaries():
    frames, raw_paths = {}, {}
    for tag, pattern in (("rt", "*RNAseq*Tumor.txt"), ("rn", "*RNAseq*Normal.txt"),
                         ("pt", "*proteomics*Tumor.txt"), ("pn", "*proteomics*Normal.txt"),
                         ("mut", "*somatic_mutation*gene_level*")):
        paths = list(RAW.glob(pattern))
        if len(paths) != 1:
            raise RuntimeError(f"Expected one source file for {tag}, found {len(paths)}")
        raw_paths[tag] = paths[0]
        d = pd.read_csv(paths[0], sep="\t")
        # A single concat avoids fragmenting the existing large matrices.
        d = d.assign(ensembl=d["idx"].str.split(".").str[0])
        frames[tag] = d
    summaries, pairs = {}, {}
    for tag, tumor_tag, normal_tag in (("rna", "rt", "rn"), ("prot", "pt", "pn")):
        t, n = frames[tumor_tag], frames[normal_tag]
        sample_ids = sorted((set(t) & set(n)) - {"idx", "ensembl"})
        pairs[tag] = sample_ids
        t = t[["ensembl"] + sample_ids].groupby("ensembl").mean()
        n = n[["ensembl"] + sample_ids].groupby("ensembl").mean()
        common = t.index.intersection(n.index)
        t, n = t.loc[common], n.loc[common]
        diff = t - n
        summaries[tag] = pd.DataFrame({
            "ensembl": common,
            f"{tag}_lfc": diff.mean(axis=1).to_numpy(),
            f"{tag}_p": stats.ttest_rel(t.to_numpy(), n.to_numpy(), axis=1,
                                        nan_policy="omit").pvalue,
            f"{tag}_n_pairs": diff.notna().sum(axis=1).to_numpy(),
        })
        # FDR applies to all measured genes tested in this layer, not just
        # the final common universe, so selection does not shrink the test family.
        summaries[tag][f"{tag}_q_bh"] = bh_qvalues(summaries[tag][f"{tag}_p"])
    mutation = frames["mut"]
    sample_ids = [c for c in mutation if c not in {"idx", "ensembl"}]
    binary = mutation[sample_ids].to_numpy()
    assert np.isin(binary[np.isfinite(binary)], [0, 1]).all()
    mut = mutation[["ensembl"] + sample_ids].groupby("ensembl").mean()
    summaries["genomic"] = pd.DataFrame({"ensembl": mut.index,
                                         "mutation_freq": mut.mean(axis=1).to_numpy()})
    source_info = {tag: {"path": p.relative_to(ROOT).as_posix(), "sha256": sha256(p),
                         "rows": len(frames[tag]), "samples": len(frames[tag].columns) - 2}
                   for tag, p in raw_paths.items()}
    return summaries, pairs, source_info


def ensure_mapping(common_ids):
    if not MAPPING.exists():
        if RESPONSE.exists():
            saved = json.loads(RESPONSE.read_text(encoding="utf-8"))
        else:
            responses = []
            for start in range(0, len(common_ids), 1000):
                subset = common_ids[start:start + 1000]
                response = requests.post("https://mygene.info/v3/query", data={
                    "q": ",".join(subset), "scopes": "ensembl.gene",
                    "fields": "symbol", "species": "human"}, timeout=(15, 90))
                response.raise_for_status()
                value = response.json()
                if not isinstance(value, list):
                    raise RuntimeError("Unexpected MyGene annotation response")
                responses.extend(value)
                print(f"MyGene: annotated queries {start + 1}-{start + len(subset)}", flush=True)
            saved = {"endpoint": "https://mygene.info/v3/query", "scopes": "ensembl.gene",
                     "fields": "symbol", "species": "human", "queried_ids": common_ids,
                     "response": responses}
            with RESPONSE.open("x", encoding="utf-8") as handle:
                json.dump(saved, handle, indent=2)
                handle.write("\n")
        assert saved["queried_ids"] == common_ids
        mapping = pd.DataFrame([{"ensembl": x.get("query"), "gene": x.get("symbol")}
                                for x in saved["response"] if x.get("symbol") is not None])
        mapping = mapping.drop_duplicates("ensembl")
        mapping.to_csv(MAPPING, sep="\t", index=False, lineterminator="\n")
    mapping = pd.read_csv(MAPPING, sep="\t")
    # Annotation services can change. Retain genuine fallback responses rather
    # than deleting a formerly mapped gene and silently changing every rank.
    saved = json.loads(RESPONSE.read_text(encoding="utf-8"))
    unresolved = [x["query"] for x in saved["response"] if x.get("symbol") is None]
    if FALLBACK_RESPONSE.exists():
        fallbacks = json.loads(FALLBACK_RESPONSE.read_text(encoding="utf-8"))
    else:
        fallbacks = []
        for identifier in unresolved:
            url = "https://rest.ensembl.org/lookup/id/" + identifier
            response = requests.get(url, params={"content-type": "application/json"}, timeout=(15, 45))
            try:
                payload = response.json()
            except ValueError:
                payload = {"raw_response_text": response.text}
            fallbacks.append({"identifier": identifier, "url": response.url,
                              "status_code": response.status_code, "response": payload})
        FALLBACK_RESPONSE.write_text(json.dumps(fallbacks, indent=2) + "\n", encoding="utf-8")
    rescued = [{"ensembl": x["identifier"], "gene": x["response"]["display_name"]}
               for x in fallbacks if x["status_code"] == 200 and x["response"].get("display_name")]
    rescue = pd.DataFrame(rescued, columns=["ensembl", "gene"])
    if not rescue.empty:
        mapping = pd.concat([mapping, rescue], ignore_index=True).drop_duplicates("ensembl")
        # Restore queried-ID order, matching the original annotation procedure.
        mapping = mapping.set_index("ensembl").reindex([i for i in common_ids if i in set(mapping.ensembl)]).reset_index()
        mapping.to_csv(MAPPING, sep="\t", index=False, lineterminator="\n")
    assert mapping["ensembl"].is_unique and mapping["gene"].notna().all()
    assert set(mapping["ensembl"]) <= set(common_ids)
    return mapping


def recover_integrated():
    summaries, pairs, sources = raw_summaries()
    ids = sorted(set(summaries["rna"]["ensembl"]) & set(summaries["prot"]["ensembl"])
                 & set(summaries["genomic"]["ensembl"]))
    mapping = ensure_mapping(ids)
    # Reproduce the professor assignment's symbol-level merges exactly.
    tx = summaries["rna"].merge(mapping, on="ensembl")
    pr = summaries["prot"].merge(mapping, on="ensembl")
    gw = summaries["genomic"].merge(mapping, on="ensembl")
    df = tx[["gene", "rna_lfc", "rna_p", "rna_n_pairs", "rna_q_bh", "ensembl"]].merge(
        pr[["gene", "prot_lfc", "prot_p", "prot_n_pairs", "prot_q_bh"]], on="gene").merge(
        gw[["gene", "mutation_freq"]], on="gene")
    assert df["gene"].is_unique, "Symbol collisions would change the original gene-level join"
    df["concordant"] = np.sign(df["rna_lfc"]) == np.sign(df["prot_lfc"])
    df["transcriptomic"], df["proteomic"], df["genomic"] = (
        df["rna_lfc"].abs(), df["prot_lfc"].abs(), df["mutation_freq"])
    for column, destination in (("transcriptomic", "tx_rank"), ("proteomic", "prot_rank"),
                                ("genomic", "gen_rank")):
        df[destination] = df[column].rank(method="average", pct=True)
    # Match the original helper's summation order, rather than changing scoring.
    df["score"] = sum((1 / 3) * df[c].rank(method="average", pct=True)
                      for c in ("transcriptomic", "proteomic", "genomic"))
    ranked = df.sort_values("score", ascending=False).reset_index(drop=True)
    ranked["module2_rank"] = np.arange(1, len(ranked) + 1)
    ranked["evidence_rna"] = ((ranked["rna_q_bh"] <= .05) &
                               (ranked["rna_lfc"].abs() >= 1) & (ranked["rna_n_pairs"] >= 10))
    ranked["evidence_protein"] = ((ranked["prot_q_bh"] <= .05) &
                                   (ranked["prot_lfc"].abs() >= 1) & (ranked["prot_n_pairs"] >= 10))
    ranked["evidence_genomic"] = ranked["mutation_freq"] >= .05
    original = pd.read_csv(ROOT / "targets_luad.csv")
    computed = ranked.head(15)[original.columns]
    pd.testing.assert_frame_equal(original.reset_index(drop=True), computed.reset_index(drop=True),
                                  check_exact=False, rtol=1e-10, atol=1e-12)
    assert len(ranked) == 6758 and int(ranked["concordant"].sum()) == 4918
    return ranked, {"raw_sources": sources, "common_ensembl_before_mapping": len(ids),
                    "mapped_ensembl": len(mapping), "rna_pair_count": len(pairs["rna"]),
                    "protein_pair_count": len(pairs["prot"]),
                    "mapping_sha256": sha256(MAPPING), "mapping_response_sha256": sha256(RESPONSE),
                    "ensembl_fallback_response_sha256": sha256(FALLBACK_RESPONSE)}


def make_integration_figure(ranked):
    top = ranked.head(15).iloc[::-1]
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 8.2), sharey=True,
                             gridspec_kw={"width_ratios": [1.5, 1, 1.4]})
    y = np.arange(len(top))
    axes[0].hlines(y, top["rna_lfc"], top["prot_lfc"], color="#C5CAD0", lw=2)
    axes[0].scatter(top["rna_lfc"], y, color="#2166AC", s=42, label="RNA")
    axes[0].scatter(top["prot_lfc"], y, color="#B2182B", s=42, label="Protein")
    axes[0].axvline(0, color="0.5", lw=.8)
    labels = [f"{int(r.module2_rank):2d}. {r.gene}" +
              ("  [discordant]" if not r.concordant else "") for r in top.itertuples()]
    axes[0].set_yticks(y, labels)
    axes[0].set_xlabel("Mean paired tumor-minus-normal (log2)")
    axes[0].set_title("A. Expression direction and magnitude", loc="left", fontweight="bold")
    axes[0].legend(frameon=False, loc="lower right")
    axes[1].barh(y, top["mutation_freq"] * 100, color="#4D8C70", height=.65)
    axes[1].set_xlabel("Tumors with mutation (%)\n108 mutation-profiled tumors")
    axes[1].set_title("B. Genomic evidence", loc="left", fontweight="bold")
    for i, value in enumerate(top["mutation_freq"] * 100):
        axes[1].text(value + .25, i, f"{value:.1f}", va="center", fontsize=9)
    axes[1].set_xlim(0, 25)
    left = np.zeros(len(top))
    for col, color, name in (("tx_rank", "#2166AC", "RNA / 3"),
                              ("prot_rank", "#B2182B", "Protein / 3"),
                              ("gen_rank", "#4D8C70", "Mutation / 3")):
        width = top[col].to_numpy() / 3
        axes[2].barh(y, width, left=left, height=.65, color=color, label=name)
        left += width
    axes[2].set_xlim(0, 1.05)
    axes[2].set_xlabel("Equal-weight rank-percentile score")
    axes[2].set_title("C. Original integrated prioritization", loc="left", fontweight="bold")
    axes[2].legend(frameon=False, loc="lower center", bbox_to_anchor=(.5, 1.06),
                   ncol=3, fontsize=8)
    fig.suptitle("LUAD multi-omics integration: original top 15 candidates", fontsize=16, y=.98)
    fig.text(.02, .035, "Ranks are preserved from targets_luad.csv. Magnitudes nominate candidates; they do not establish therapeutic benefit.\n"
             "Protein coverage varies by gene (see integrated table); ZNF536's protein signal uses 15 valid paired measurements.",
             fontsize=10, ha="left")
    fig.tight_layout(rect=(0, .105, 1, .95), w_pad=2)
    fig.savefig(FIGURES / "luad_multiomics_integration.png", dpi=220)
    fig.savefig(FIGURES / "luad_multiomics_integration.pdf")
    plt.close(fig)


def make_venn(ranked):
    sets = {"RNA": set(ranked.loc[ranked.evidence_rna, "gene"]),
            "Protein": set(ranked.loc[ranked.evidence_protein, "gene"]),
            "Genomic": set(ranked.loc[ranked.evidence_genomic, "gene"])}
    a, b, c = sets.values()
    regions = {"100": a - b - c, "010": b - a - c, "001": c - a - b,
               "110": (a & b) - c, "101": (a & c) - b,
               "011": (b & c) - a, "111": a & b & c}
    outside = set(ranked.gene) - (a | b | c)
    fig, ax = plt.subplots(figsize=(9.2, 8.4))
    colors = ("#2166AC", "#B2182B", "#4D8C70")
    for center, color in zip(((.38, .60), (.62, .60), (.50, .38)), colors):
        ax.add_patch(Circle(center, .27, facecolor=color, edgecolor=color, alpha=.24, lw=2))
    positions = {"100": (.25, .64), "010": (.75, .64), "001": (.50, .16),
                 "110": (.50, .77), "101": (.32, .39), "011": (.68, .39), "111": (.50, .53)}
    for key, values in regions.items():
        ax.text(*positions[key], f"{len(values):,}", ha="center", va="center", fontsize=17, fontweight="bold")
    ax.text(.18, .92, f"RNA\nn = {len(a):,}", ha="center", va="center", fontsize=13, color=colors[0])
    ax.text(.82, .92, f"Protein\nn = {len(b):,}", ha="center", va="center", fontsize=13, color=colors[1])
    ax.text(.50, .015, f"Genomic recurrence\nn = {len(c):,}", ha="center", va="center", fontsize=13, color=colors[2])
    ax.set(xlim=(.02, .98), ylim=(-.04, 1.02), aspect="equal")
    ax.axis("off")
    fig.suptitle("LUAD: overlap of thresholded positive evidence", fontsize=16, y=.98)
    fig.text(.08, .105, "RNA / protein: BH q ≤ 0.05, |mean log2 difference| ≥ 1, and ≥ 10 valid pairs.\n"
             "Genomic: mutation in ≥ 5% of 108 tumors (descriptive recurrence, not a significance test).\n"
             f"Universe: 6,758 jointly measured mapped genes; outside all sets: {len(outside):,}.\n"
             "Region counts are exact; circles are schematic, not area-proportional. Thresholds do not alter target ranks.",
             fontsize=10, ha="left")
    fig.tight_layout(rect=(0, .19, 1, .96))
    fig.savefig(FIGURES / "luad_evidence_venn.png", dpi=220)
    fig.savefig(FIGURES / "luad_evidence_venn.pdf")
    plt.close(fig)
    membership = ranked[["gene", "ensembl", "evidence_rna", "evidence_protein", "evidence_genomic"]].copy()
    membership["region_RNA_Protein_Genomic"] = membership[["evidence_rna", "evidence_protein", "evidence_genomic"]].astype(int).astype(str).agg("".join, axis=1)
    membership.to_csv(DATA / "luad_evidence_set_membership.csv", index=False)
    return {"set_sizes": {key: len(value) for key, value in sets.items()},
            "exclusive_regions_RNA_Protein_Genomic": {key: len(value) for key, value in regions.items()},
            "outside_all_sets": len(outside)}


def main():
    DATA.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    # Persist threshold specification before creating any memberships/overlap.
    threshold_path = DATA / "luad_evidence_thresholds.json"
    if threshold_path.exists():
        assert json.loads(threshold_path.read_text()) == EVIDENCE_THRESHOLDS
    else:
        threshold_path.write_text(json.dumps(EVIDENCE_THRESHOLDS, indent=2) + "\n", encoding="utf-8")
    original_target_hash = sha256(ROOT / "targets_luad.csv")
    ranked, provenance = recover_integrated()
    ranked.to_csv(DATA / "luad_integrated_gene_level.csv", index=False)
    known = ranked[ranked.gene.isin(["TP53", "KRAS", "EGFR", "STK11", "KEAP1", "BRAF", "MET", "ERBB2", "ALK", "NF1", "ZNF536"])]
    known.to_csv(DATA / "luad_known_and_discordant_genes.csv", index=False)
    make_integration_figure(ranked)
    overlap = make_venn(ranked)
    assert sha256(ROOT / "targets_luad.csv") == original_target_hash
    summary = {"joined_genes": len(ranked), "concordant": int(ranked.concordant.sum()),
               "discordant": int((~ranked.concordant).sum()), "exported_original_targets": 15,
               "original_targets_sha256_unchanged": original_target_hash,
               "scoring": "Equal mean of rank-percentiles of absolute RNA and protein effects and mutation frequency; no QC re-ranking",
               "evidence_thresholds": EVIDENCE_THRESHOLDS, **overlap, **provenance,
               "numpy_version": np.__version__, "pandas_version": pd.__version__,
               "gene_details": known.where(pd.notna(known), None).to_dict(orient="records")}
    (DATA / "luad_portfolio_verification.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("joined_genes", "concordant", "discordant", "set_sizes", "exclusive_regions_RNA_Protein_Genomic", "outside_all_sets")}, indent=2))
    print(known[["gene", "module2_rank", "rna_lfc", "prot_lfc", "mutation_freq", "prot_n_pairs", "prot_q_bh"]].to_string(index=False))


if __name__ == "__main__":
    main()
