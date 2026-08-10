from __future__ import annotations

import itertools
import json
import math
from pathlib import Path
from typing import Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon


BASE = Path(__file__).resolve().parents[1]
PROCESSED = BASE / "processed"
RESULTS = BASE / "results"
FIGURES = BASE / "figures"


def node_map(edges: pd.DataFrame) -> dict[str, set[str]]:
    return edges.groupby("canonical_herb_id")["node_id"].apply(set).to_dict()


def union_nodes(herbs: list[str] | tuple[str, ...], mapping: Mapping[str, set[str]]) -> set[str]:
    out: set[str] = set()
    for herb in herbs:
        out.update(mapping.get(herb, set()))
    return out


def exact_curve(herbs: list[str], mapping: Mapping[str, set[str]]) -> pd.DataFrame:
    n = len(herbs)
    full = union_nodes(herbs, mapping)
    support = [sum(node in mapping.get(h, set()) for h in herbs) for node in full]
    rows = []
    for removed in range(n + 1):
        den = math.comb(n, removed)
        retention = [1.0 - (math.comb(n - r, removed) / den if removed <= n - r else 0.0) for r in support]
        rows.append({"removed_n": removed, "fraction_removed": removed / n, "mean_retention": float(np.mean(retention))})
    return pd.DataFrame(rows)


def fdi(curve: pd.DataFrame) -> float:
    auc = np.trapezoid(curve["mean_retention"], curve["fraction_removed"])
    return 1.0 - float(auc)


def shapley_union(herbs: list[str], mapping: Mapping[str, set[str]]) -> pd.DataFrame:
    full = union_nodes(herbs, mapping)
    values = {h: 0.0 for h in herbs}
    unique = {h: 0.0 for h in herbs}
    redundant = {h: 0.0 for h in herbs}
    for node in full:
        supporters = [h for h in herbs if node in mapping.get(h, set())]
        increment = 1.0 / (len(supporters) * len(full))
        for herb in supporters:
            values[herb] += increment
            (unique if len(supporters) == 1 else redundant)[herb] += increment
    out = pd.DataFrame({"herb": herbs, "shapley": [values[h] for h in herbs], "unique_component": [unique[h] for h in herbs], "redundant_component": [redundant[h] for h in herbs]})
    assert np.isclose(out["shapley"].sum(), 1.0)
    return out.sort_values("shapley", ascending=False)


def gini(values: np.ndarray) -> float:
    values = np.sort(np.asarray(values, dtype=float))
    if values.sum() == 0:
        return 0.0
    n = len(values)
    return float((2 * np.sum((np.arange(1, n + 1)) * values) / (n * values.sum())) - (n + 1) / n)


def topology_metrics(herbs: list[str], mapping: Mapping[str, set[str]], shapley: pd.DataFrame) -> dict[str, float]:
    full = union_nodes(herbs, mapping)
    counts = np.array([sum(node in mapping.get(h, set()) for h in herbs) for node in full], dtype=float)
    p = counts / counts.sum()
    entropy = -float(np.sum(p * np.log(p)))
    loho = [len(union_nodes([x for x in herbs if x != h], mapping)) / len(full) for h in herbs]
    return {
        "n_nodes": len(full),
        "mean_herb_support": float(counts.mean()),
        "multiherb_fraction": float(np.mean(counts >= 2)),
        "redundancy_ratio": float(min(loho)),
        "mean_loho_retention": float(np.mean(loho)),
        "normalized_entropy": entropy / math.log(len(full)) if len(full) > 1 else 0.0,
        "entropy_equivalent_nodes": float(np.exp(entropy)),
        "shapley_hhi": float(np.square(shapley["shapley"]).sum()),
        "shapley_gini": gini(shapley["shapley"].to_numpy()),
        "max_shapley": float(shapley["shapley"].max()),
    }


def enumerate_subsets(formula_id: str, node_type: str, herbs: list[str], mapping: Mapping[str, set[str]]) -> list[dict]:
    full = union_nodes(herbs, mapping)
    rows = []
    for retained_n in range(1, len(herbs) + 1):
        for subset in itertools.combinations(herbs, retained_n):
            nodes = union_nodes(subset, mapping)
            rows.append({"formula_id": formula_id, "node_type": node_type, "formula_size": len(herbs), "retained_n": retained_n, "removed_n": len(herbs) - retained_n, "subset": "|".join(subset), "node_count": len(nodes), "retention": len(nodes) / len(full)})
    return rows


def robust_model(summary: pd.DataFrame, outcome: str) -> tuple[pd.DataFrame, str]:
    frame = summary[summary["node_type"] == outcome].copy()
    x = np.column_stack(
        [
            np.ones(len(frame)),
            frame["formula_size"].to_numpy(dtype=float),
            np.log1p(frame["n_nodes"].to_numpy(dtype=float)),
            frame["multiherb_fraction"].to_numpy(dtype=float),
        ]
    )
    y = frame["fdi"].to_numpy(dtype=float)
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residual = y - x @ beta
    leverage = np.sum((x @ xtx_inv) * x, axis=1)
    scaled = residual / np.clip(1.0 - leverage, 1e-8, None)
    meat = x.T @ (x * np.square(scaled)[:, None])
    covariance = xtx_inv @ meat @ xtx_inv
    se = np.sqrt(np.clip(np.diag(covariance), 0, None))
    from scipy.stats import norm as normal_distribution

    z = np.divide(beta, se, out=np.full_like(beta, np.nan), where=se > 0)
    p = 2 * normal_distribution.sf(np.abs(z))
    terms = ["Intercept", "formula_size", "log1p_n_nodes", "multiherb_fraction"]
    table = pd.DataFrame(
        {
            "term": terms,
            "estimate": beta,
            "se_hc3": se,
            "p_value": p,
            "ci_low": beta - 1.96 * se,
            "ci_high": beta + 1.96 * se,
        }
    )
    fitted = x @ beta
    r2 = 1.0 - np.square(y - fitted).sum() / np.square(y - y.mean()).sum()
    text = f"Manual OLS with HC3 covariance; n={len(y)}; R2={r2:.6f}\n" + table.to_string(index=False)
    return table, text


def make_figure(summary_wide: pd.DataFrame, identity: pd.DataFrame, gates: pd.DataFrame) -> None:
    colors = {"compound": "#E76F51", "target": "#277DA1"}
    fig = plt.figure(figsize=(14, 10), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.9, 1.1])

    ax = fig.add_subplot(gs[0, 0])
    candidate_n = int((identity["formula_id"] != "F000").sum())
    included_n = int(((identity["formula_id"] != "F000") & identity["included_identity_gate"].astype(bool)).sum())
    ax.axis("off")
    boxes = [(0.05, 0.62, f"Candidate formulae\n{candidate_n}"), (0.38, 0.62, f"Identity-qualified\n{included_n}"), (0.71, 0.62, f"Topology benchmark\n{included_n} + SQRCD")]
    for x, y, label in boxes:
        ax.add_patch(plt.Rectangle((x, y), 0.24, 0.22, facecolor="#E8F1F8", edgecolor="#28536B", lw=1.5, transform=ax.transAxes))
        ax.text(x + 0.12, y + 0.11, label, ha="center", va="center", fontsize=10, transform=ax.transAxes)
    ax.annotate("", xy=(0.38, 0.73), xytext=(0.29, 0.73), arrowprops=dict(arrowstyle="->", color="#555"), xycoords=ax.transAxes)
    ax.annotate("", xy=(0.71, 0.73), xytext=(0.62, 0.73), arrowprops=dict(arrowstyle="->", color="#555"), xycoords=ax.transAxes)
    ax.text(0.5, 0.28, f"Excluded at identity gate: {candidate_n - included_n}\nPrespecified Jaccard < 0.80 or size outside 4–25", ha="center", va="center", fontsize=9, color="#7A3E3E", transform=ax.transAxes)
    ax.set_title("A  Cross-formula inclusion", loc="left", fontweight="bold")

    ax = fig.add_subplot(gs[0, 1])
    order = summary_wide.sort_values("target_fdi")["formula_name_en"].tolist()
    y = np.arange(len(order))
    temp = summary_wide.set_index("formula_name_en").loc[order]
    ax.hlines(y, temp["target_fdi"], temp["compound_fdi"], color="#B7C4CE", lw=2)
    ax.scatter(temp["compound_fdi"], y, color=colors["compound"], label="Compound FDI", s=35)
    ax.scatter(temp["target_fdi"], y, color=colors["target"], label="Target FDI", s=35)
    if "SQRCD" in order:
        sq = order.index("SQRCD")
        ax.axhspan(sq - 0.45, sq + 0.45, color="#F4D35E", alpha=0.25)
    ax.set_yticks(y, order, fontsize=8)
    ax.set_xlim(0, max(0.7, float(temp[["compound_fdi", "target_fdi"]].max().max()) + 0.05))
    ax.set_xlabel("Formula Discrimination Index")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("B  Harmonized FDI across formulae", loc="left", fontweight="bold")

    ax = fig.add_subplot(gs[1, 0])
    scatter = ax.scatter(summary_wide["target_multiherb_fraction"], summary_wide["target_fdi"], s=25 + summary_wide["formula_size"] * 8, c=summary_wide["formula_size"], cmap="viridis", alpha=0.85, edgecolor="white", linewidth=0.7)
    sq = summary_wide[summary_wide["formula_id"] == "F000"]
    if not sq.empty:
        ax.scatter(sq["target_multiherb_fraction"], sq["target_fdi"], marker="*", s=260, c="#D62828", edgecolor="black", linewidth=0.8, label="SQRCD")
        ax.legend(frameon=False)
    ax.set_xlabel("Multi-herb target fraction")
    ax.set_ylabel("Target FDI")
    ax.set_title("C  Redundancy and formula discrimination", loc="left", fontweight="bold")
    fig.colorbar(scatter, ax=ax, label="Formula size")

    ax = fig.add_subplot(gs[1, 1])
    state_code = {"qualified": 3, "partially qualified": 2, "not estimable": 1, "stopped": 0}
    matrix = gates.set_index("formula_name_en").loc[order]
    arr = matrix.drop(columns=["formula_id"]).replace(state_code).to_numpy(dtype=float)
    cmap = plt.matplotlib.colors.ListedColormap(["#C44E52", "#B8B8B8", "#F2C14E", "#4C9F70"])
    ax.imshow(arr, aspect="auto", cmap=cmap, vmin=0, vmax=3)
    ax.set_yticks(np.arange(len(order)), order, fontsize=8)
    ax.set_xticks(np.arange(arr.shape[1]), [x.replace("_", "\n") for x in matrix.drop(columns=["formula_id"]).columns], fontsize=7)
    ax.set_title("D  Gate profiles delimit benchmark scope", loc="left", fontweight="bold")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=cmap(i), label=s) for s, i in state_code.items()], loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=7)

    fig.suptitle("External benchmarking of EQF topology modules", fontsize=16, fontweight="bold")
    fig.savefig(FIGURES / "Figure6_external_benchmark.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "Figure6_external_benchmark.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    formulas = pd.read_csv(PROCESSED / "formulas_locked.csv")
    mapping_audit = pd.read_csv(PROCESSED / "herb_mapping_audit.csv")
    identity = pd.read_csv(PROCESSED / "formula_identity_audit.csv")
    compound_map = node_map(pd.read_csv(PROCESSED / "herb_compound_edges.csv"))
    target_map = node_map(pd.read_csv(PROCESSED / "herb_target_edges.csv"))

    summaries, curves, shapleys, subsets, loho_rows = [], [], [], [], []
    for formula_id, group in formulas.groupby("formula_id"):
        herbs = sorted(group["canonical_herb_id"].unique())
        name = group.iloc[0]["formula_name_en"]
        mapping_rate = float(mapping_audit.loc[mapping_audit["formula_id"] == formula_id, "mapping_rate"].iloc[0])
        for node_type, mapping in [("compound", compound_map), ("target", target_map)]:
            curve = exact_curve(herbs, mapping)
            curve.insert(0, "formula_id", formula_id)
            curve.insert(1, "node_type", node_type)
            curves.append(curve)
            shapley = shapley_union(herbs, mapping)
            shapley.insert(0, "formula_id", formula_id)
            shapley.insert(1, "node_type", node_type)
            shapleys.append(shapley)
            metrics = topology_metrics(herbs, mapping, shapley)
            summaries.append({"formula_id": formula_id, "formula_name_en": name, "node_type": node_type, "formula_size": len(herbs), "mapping_rate": mapping_rate, "fdi": fdi(curve), **metrics})
            full_n = len(union_nodes(herbs, mapping))
            for herb in herbs:
                retained = len(union_nodes([h for h in herbs if h != herb], mapping)) / full_n
                loho_rows.append({"formula_id": formula_id, "node_type": node_type, "removed_herb": herb, "retention": retained, "loss": 1 - retained})
            subsets.extend(enumerate_subsets(formula_id, node_type, herbs, mapping))

    summary = pd.DataFrame(summaries)
    curve_df = pd.concat(curves, ignore_index=True)
    shapley_df = pd.concat(shapleys, ignore_index=True)
    subset_df = pd.DataFrame(subsets)
    loho_df = pd.DataFrame(loho_rows)
    summary.to_csv(RESULTS / "external_benchmark_summary_long.csv", index=False)
    curve_df.to_csv(RESULTS / "formula_retention_curves.csv", index=False)
    shapley_df.to_csv(RESULTS / "formula_shapley.csv", index=False)
    subset_df.to_csv(RESULTS / "complete_subformula_lattice.csv.gz", index=False, compression="gzip")
    loho_df.to_csv(RESULTS / "formula_loho.csv", index=False)

    wide = summary.pivot(index=["formula_id", "formula_name_en", "formula_size", "mapping_rate"], columns="node_type", values=["fdi", "n_nodes", "multiherb_fraction", "mean_loho_retention", "max_shapley", "shapley_hhi", "normalized_entropy"]).reset_index()
    wide.columns = ["_".join([str(x) for x in col if str(x)]) if isinstance(col, tuple) else col for col in wide.columns]
    rename = {}
    for col in wide.columns:
        for metric in ["fdi", "n_nodes", "multiherb_fraction", "mean_loho_retention", "max_shapley", "shapley_hhi", "normalized_entropy"]:
            if col.startswith(metric + "_"):
                rename[col] = col[len(metric) + 1 :] + "_" + metric
    wide = wide.rename(columns=rename)
    wide["size_group"] = pd.cut(wide["formula_size"], bins=[0, 6, 10, 15, 25], labels=["4-6", "7-10", "11-15", "16-25"], include_lowest=True)
    wide["comparison_band"] = pd.cut(wide["formula_size"], bins=[0, 6, 10, 20, 25], labels=["4-6", "7-10", "11-20", "21-25"], include_lowest=True)
    wide["target_fdi_percentile_within_band"] = wide.groupby("comparison_band", observed=True)["target_fdi"].rank(pct=True)
    wide.to_csv(RESULTS / "external_benchmark_summary.csv", index=False)

    compound = wide["compound_fdi"].to_numpy()
    target = wide["target_fdi"].to_numpy()
    paired = wilcoxon(compound, target, alternative="two-sided")
    rho, rho_p = spearmanr(wide["target_multiherb_fraction"], wide["target_fdi"])
    model_tables = []
    model_text = []
    for node_type in ["compound", "target"]:
        table, text = robust_model(summary, node_type)
        table.insert(0, "node_type", node_type)
        model_tables.append(table)
        model_text.append(f"=== {node_type.upper()} FDI MODEL ===\n{text}")
    pd.concat(model_tables, ignore_index=True).to_csv(RESULTS / "fdi_hc3_models.csv", index=False)
    (RESULTS / "fdi_hc3_models.txt").write_text("\n\n".join(model_text), encoding="utf-8")

    gate_rows = []
    for row in wide.itertuples():
        gate_rows.append({"formula_id": row.formula_id, "formula_name_en": row.formula_name_en, "identity": "qualified", "chemical_topology": "partially qualified", "target_topology": "partially qualified", "preparation_chemistry": "partially qualified" if row.formula_id == "F000" else "not estimable", "human_context": "not estimable", "genetics": "not estimable", "perturbation": "not estimable", "structural_eligibility": "stopped"})
    gates = pd.DataFrame(gate_rows)
    gates.to_csv(RESULTS / "benchmark_gate_profiles.csv", index=False)
    make_figure(wide, identity, gates)

    sq = wide[wide["formula_id"] == "F000"].iloc[0]
    report = {
        "included_formulae_total": int(len(wide)),
        "external_formulae": int(len(wide) - 1),
        "formula_size_range": [int(wide["formula_size"].min()), int(wide["formula_size"].max())],
        "median_mapping_rate": float(wide["mapping_rate"].median()),
        "compound_fdi_range": [float(wide["compound_fdi"].min()), float(wide["compound_fdi"].max())],
        "target_fdi_range": [float(wide["target_fdi"].min()), float(wide["target_fdi"].max())],
        "target_lower_than_compound_n": int((wide["target_fdi"] < wide["compound_fdi"]).sum()),
        "paired_wilcoxon_statistic": float(paired.statistic),
        "paired_wilcoxon_p": float(paired.pvalue),
        "target_fdi_vs_redundancy_spearman_rho": float(rho),
        "target_fdi_vs_redundancy_p": float(rho_p),
        "sqrcd_harmonized_compound_fdi": float(sq["compound_fdi"]),
        "sqrcd_harmonized_target_fdi": float(sq["target_fdi"]),
        "sqrcd_target_fdi_percentile_11_20": float(sq["target_fdi_percentile_within_band"]),
        "complete_subformula_rows_both_node_types": int(len(subset_df)),
        "interpretive_boundary": "topology modules externally benchmarked; complete eight-gate EQF not externally validated",
    }
    (RESULTS / "benchmark_key_results.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
