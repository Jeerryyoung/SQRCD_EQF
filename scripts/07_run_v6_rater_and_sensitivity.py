from __future__ import annotations

import json
import math
import re
import unicodedata
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import cohen_kappa_score


ROOT = Path(__file__).resolve().parents[1]
DATE = ROOT / "Date"
EXT = ROOT / "external_benchmark"
SUPP = ROOT / "06_supplementary"
EQF = ROOT / "02_eqf"
ROB = ROOT / "02_robustness"
FIG = ROOT / "05_figures"
for folder in (SUPP, EQF, ROB, FIG):
    folder.mkdir(exist_ok=True)


def find_one(pattern: str) -> Path:
    hits = sorted(DATE.glob(pattern))
    if len(hits) != 1:
        raise RuntimeError(f"Expected one file for {pattern!r}; found {hits}")
    return hits[0]


blank_path = find_one("1-*.xlsx")
rater_a_path = find_one("2-*.xlsx")
rater_b_path = find_one("*ASSESSED.xlsx")


def read_rater(path: Path) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name="Rater_Form", dtype=object)


def read_blank(path: Path) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name="Second_Rater_Assessment", header=2, dtype=object)


blank = read_blank(blank_path)
a = read_rater(rater_a_path)
b = read_rater(rater_b_path)
key = ["dossier_id", "rubric_code"]


def audit_row(label: str, frame: pd.DataFrame, rating_col: str, initials_col: str, date_col: str) -> dict:
    required = 750
    valid = {"qualified", "partially qualified", "not estimable", "stopped"}
    states = frame[rating_col].dropna().astype(str).str.strip().str.lower()
    return {
        "file_role": label,
        "rows": len(frame),
        "formulae": frame[frame.columns[0]].nunique(),
        "rubric_items": frame[frame.columns[3]].nunique(),
        "duplicate_formula_item_keys": int(frame.duplicated([frame.columns[0], frame.columns[3]]).sum()),
        "completed_ratings": int(frame[rating_col].notna().sum()),
        "completion_rate": float(frame[rating_col].notna().mean()),
        "invalid_state_rows": int((~states.isin(valid)).sum()),
        "distinct_initials": "|".join(sorted(frame[initials_col].dropna().astype(str).unique())),
        "distinct_dates": "|".join(sorted(frame[date_col].dropna().astype(str).unique())),
        "usable_for_agreement": bool(len(frame) == required and frame[rating_col].notna().all() and frame.duplicated([frame.columns[0], frame.columns[3]]).sum() == 0 and states.isin(valid).all()),
    }


audit = pd.DataFrame([
    audit_row("blank_template", blank, "Second Rater Rating", "Rater Initials", "Rating Date"),
    audit_row("completed_rating_A", a, "rater_state", "rater_initials", "rating_date"),
    audit_row("completed_rating_B", b, "rater_state", "rater_initials", "rating_date"),
])

merged = a.merge(b, on=key, suffixes=("_A", "_B"), validate="one_to_one")
if len(merged) != 750:
    raise AssertionError("The two completed ratings do not share all 750 keys.")

evidence_columns = ["formula_name_en", "formula_name_zh", "domain", "field", "criterion", "organizer_availability", "raw_evidence", "source_locator"]
for col in evidence_columns:
    if not merged[f"{col}_A"].fillna("").equals(merged[f"{col}_B"].fillna("")):
        raise AssertionError(f"Evidence column differs between raters: {col}")

states = ["qualified", "partially qualified", "not estimable", "stopped"]
for col in ["rater_state_A", "rater_state_B"]:
    merged[col] = merged[col].astype(str).str.strip().str.lower()
    if not set(merged[col]).issubset(states):
        raise AssertionError(f"Invalid states in {col}: {set(merged[col]) - set(states)}")


def collapse_unavailable(series: pd.Series) -> pd.Series:
    return series.replace({"stopped": "unavailable", "not estimable": "unavailable"})


def metric_rows(label: str, x: pd.Series, y: pd.Series) -> list[dict]:
    labels = sorted(set(x) | set(y))
    return [
        {"analysis": label, "metric": "n", "value": len(x)},
        {"analysis": label, "metric": "agreement_n", "value": int((x == y).sum())},
        {"analysis": label, "metric": "raw_agreement", "value": float((x == y).mean())},
        {"analysis": label, "metric": "cohen_kappa_unweighted", "value": float(cohen_kappa_score(x, y, labels=labels))},
    ]


agreement_rows = metric_rows("four_state", merged.rater_state_A, merged.rater_state_B)
collapsed_a = collapse_unavailable(merged.rater_state_A)
collapsed_b = collapse_unavailable(merged.rater_state_B)
agreement_rows += metric_rows("collapsed_unavailable", collapsed_a, collapsed_b)

# Linear weighting is a prespecified sensitivity only. Not-estimable and stopped are process states,
# so unweighted kappa remains the primary four-state statistic.
ordinal = {"stopped": 0, "not estimable": 1, "partially qualified": 2, "qualified": 3}
agreement_rows.append({
    "analysis": "four_state_ordered_sensitivity",
    "metric": "cohen_kappa_linear_weighted",
    "value": float(cohen_kappa_score(merged.rater_state_A.map(ordinal), merged.rater_state_B.map(ordinal), weights="linear")),
})
agreement_summary = pd.DataFrame(agreement_rows)


domain_rows = []
for domain, g in merged.groupby("domain_A", sort=False):
    xa, xb = g.rater_state_A, g.rater_state_B
    for analysis, xx, yy in [
        ("four_state", xa, xb),
        ("collapsed_unavailable", collapse_unavailable(xa), collapse_unavailable(xb)),
    ]:
        labels = sorted(set(xx) | set(yy))
        kappa = float(cohen_kappa_score(xx, yy, labels=labels)) if len(labels) > 1 else np.nan
        domain_rows.append({
            "domain": domain,
            "analysis": analysis,
            "n": len(g),
            "agreement_n": int((xx == yy).sum()),
            "raw_agreement": float((xx == yy).mean()),
            "cohen_kappa_unweighted": kappa,
        })
domain_agreement = pd.DataFrame(domain_rows)

formula_rows = []
for fid, g in merged.groupby("dossier_id", sort=False):
    for analysis, xx, yy in [
        ("four_state", g.rater_state_A, g.rater_state_B),
        ("collapsed_unavailable", collapse_unavailable(g.rater_state_A), collapse_unavailable(g.rater_state_B)),
    ]:
        labels = sorted(set(xx) | set(yy))
        kappa = float(cohen_kappa_score(xx, yy, labels=labels)) if len(labels) > 1 else np.nan
        formula_rows.append({
            "formula_id": fid,
            "formula_name_en": g.formula_name_en_A.iloc[0],
            "analysis": analysis,
            "n": len(g),
            "agreement_n": int((xx == yy).sum()),
            "raw_agreement": float((xx == yy).mean()),
            "cohen_kappa_unweighted": kappa,
        })
formula_agreement = pd.DataFrame(formula_rows)

disagreements = merged.loc[merged.rater_state_A.ne(merged.rater_state_B), [
    "dossier_id", "formula_name_en_A", "rubric_code", "domain_A", "field_A", "organizer_availability_A", "rater_state_A", "rater_state_B", "rater_note_A", "rater_note_B"
]].rename(columns={c: c.replace("_A", "") for c in ["formula_name_en_A", "domain_A", "field_A", "organizer_availability_A"]})


def cluster_bootstrap(frame: pd.DataFrame, reps: int = 2000, seed: int = 20260802) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    groups = {fid: g for fid, g in frame.groupby("dossier_id", sort=False)}
    ids = np.array(list(groups))
    out = []
    for mode in ["four_state", "collapsed_unavailable"]:
        vals = []
        for _ in range(reps):
            sample_ids = rng.choice(ids, size=len(ids), replace=True)
            sample = pd.concat([groups[x] for x in sample_ids], ignore_index=True)
            x, y = sample.rater_state_A, sample.rater_state_B
            if mode == "collapsed_unavailable":
                x, y = collapse_unavailable(x), collapse_unavailable(y)
            labels = sorted(set(x) | set(y))
            kappa = cohen_kappa_score(x, y, labels=labels) if len(labels) > 1 else np.nan
            vals.append(((x == y).mean(), kappa))
        arr = np.asarray(vals, float)
        for metric, col in [("raw_agreement", 0), ("cohen_kappa_unweighted", 1)]:
            finite = arr[:, col][np.isfinite(arr[:, col])]
            out.append({
                "analysis": mode,
                "metric": metric,
                "replicates": reps,
                "estimate": float(agreement_summary.loc[(agreement_summary.analysis.eq(mode)) & (agreement_summary.metric.eq(metric)), "value"].iloc[0]),
                "cluster_bootstrap_q025": float(np.quantile(finite, 0.025)),
                "cluster_bootstrap_q975": float(np.quantile(finite, 0.975)),
                "seed": seed,
            })
    return pd.DataFrame(out)


bootstrap_agreement = cluster_bootstrap(merged)

# Data quality findings are evidence-use decisions, not corrections to the original workbooks.
audit["finding"] = [
    "Template structure is complete, but all 750 rating, initials and date cells are blank; not usable as a rater record.",
    "Complete 750-item record with unique formula-item keys and allowed states.",
    "Complete 750-item record with unique formula-item keys and allowed states.",
]
audit["risk_or_boundary"] = [
    "Do not count as a third evaluator.",
    "Distinct initials/date are present, but independence cannot be established from workbook metadata alone.",
    "Distinct initials/date are present, but independence cannot be established from workbook metadata alone.",
]

audit.to_csv(SUPP / "second_rater_data_quality_audit.csv", index=False, encoding="utf-8-sig")
agreement_summary.to_csv(EQF / "second_rater_agreement_summary.csv", index=False, encoding="utf-8-sig")
domain_agreement.to_csv(EQF / "second_rater_domain_agreement.csv", index=False, encoding="utf-8-sig")
formula_agreement.to_csv(EQF / "second_rater_formula_agreement.csv", index=False, encoding="utf-8-sig")
bootstrap_agreement.to_csv(EQF / "second_rater_cluster_bootstrap.csv", index=False, encoding="utf-8-sig")
disagreements.to_csv(EQF / "second_rater_disagreements.csv", index=False, encoding="utf-8-sig")
merged.to_csv(EQF / "second_rater_field_level_joined.csv", index=False, encoding="utf-8-sig")


# Compound entity harmonization sensitivity.
formulas = pd.read_csv(EXT / "processed" / "formulas_locked.csv", encoding="utf-8-sig")
edges = pd.read_csv(EXT / "processed" / "herb_compound_edges.csv", encoding="utf-8-sig")
benchmark = pd.read_csv(EXT / "results" / "external_benchmark_summary.csv", encoding="utf-8-sig")


def normalize_name(value: str) -> str:
    value = unicodedata.normalize("NFKC", str(value)).casefold().strip()
    value = value.replace("α", "alpha").replace("β", "beta").replace("γ", "gamma")
    value = re.sub(r"[\s_]+", " ", value)
    value = re.sub(r"\s*([,;:/()\[\]-])\s*", r"\1", value)
    value = re.sub(r"[^0-9a-z+\-(),;/\[\] ]+", "", value)
    return value.strip()


edges["normalized_name"] = edges.node_label.map(normalize_name)
edges = edges.loc[edges.normalized_name.ne("")].copy()


def fdi_from_map(herbs: list[str], mapping: dict[str, set[str]]) -> tuple[float, int, float]:
    full = set().union(*(mapping.get(h, set()) for h in herbs))
    n = len(herbs)
    if not full:
        return np.nan, 0, np.nan
    support = np.array([sum(node in mapping.get(h, set()) for h in herbs) for node in full], dtype=int)
    xs, ys = [], []
    for removed in range(n + 1):
        denom = math.comb(n, removed)
        retention = []
        for r in support:
            absent = 0 if removed > n - r else math.comb(n - r, removed) / denom
            retention.append(1 - absent)
        xs.append(removed / n)
        ys.append(float(np.mean(retention)))
    return 1 - float(np.trapezoid(ys, xs)), len(full), float(np.mean(support >= 2))


maps = {
    "ETCM_native_record_id": edges.groupby("canonical_herb_id").node_id.apply(set).to_dict(),
    "normalized_compound_name": edges.groupby("canonical_herb_id").normalized_name.apply(set).to_dict(),
}
entity_rows = []
for fid, g in formulas.groupby("formula_id", sort=False):
    herbs = sorted(g.canonical_herb_id.dropna().unique())
    for representation, mapping_dict in maps.items():
        fdi, n_nodes, multi = fdi_from_map(herbs, mapping_dict)
        entity_rows.append({
            "formula_id": fid,
            "formula_name_en": g.formula_name_en.iloc[0],
            "formula_size": len(herbs),
            "entity_representation": representation,
            "status": "estimated",
            "compound_nodes": n_nodes,
            "compound_fdi": fdi,
            "multiherb_fraction": multi,
            "note": "ETCM-native identifiers retained" if representation == "ETCM_native_record_id" else "Unicode, case, whitespace and punctuation normalization of ETCM compound labels; no salt or stereochemistry collapse",
        })
    for representation in ["PubChem_CID", "parent_InChIKey"]:
        entity_rows.append({
            "formula_id": fid,
            "formula_name_en": g.formula_name_en.iloc[0],
            "formula_size": len(herbs),
            "entity_representation": representation,
            "status": "not estimable",
            "compound_nodes": np.nan,
            "compound_fdi": np.nan,
            "multiherb_fraction": np.nan,
            "note": "Uniform identifier coverage is absent from the retained ETCM external-benchmark export; no cross-formula imputation was performed.",
        })
entity = pd.DataFrame(entity_rows)
entity.to_csv(ROB / "compound_entity_harmonization_sensitivity.csv", index=False, encoding="utf-8-sig")

coverage = pd.DataFrame([
    {"representation": "ETCM_native_record_id", "covered_edges": len(edges), "total_edges": len(edges), "coverage": 1.0, "status": "available"},
    {"representation": "normalized_compound_name", "covered_edges": int(edges.normalized_name.ne("").sum()), "total_edges": len(edges), "coverage": float(edges.normalized_name.ne("").mean()), "status": "available"},
    {"representation": "PubChem_CID", "covered_edges": 0, "total_edges": len(edges), "coverage": 0.0, "status": "not estimable"},
    {"representation": "parent_InChIKey", "covered_edges": 0, "total_edges": len(edges), "coverage": 0.0, "status": "not estimable"},
])
coverage.to_csv(ROB / "compound_entity_identifier_coverage.csv", index=False, encoding="utf-8-sig")


# Cross-formula leave-one-formula-out and bootstrap sensitivity.
x = benchmark.target_multiherb_fraction.to_numpy(float)
y = benchmark.target_fdi.to_numpy(float)
rho_full, p_full = spearmanr(x, y)
lofo_rows = []
for i, row in benchmark.reset_index(drop=True).iterrows():
    keep = np.arange(len(benchmark)) != i
    rho, p = spearmanr(x[keep], y[keep])
    X = np.column_stack([
        benchmark.loc[keep, "formula_size"].to_numpy(float),
        np.log1p(benchmark.loc[keep, "target_n_nodes"].to_numpy(float)),
        benchmark.loc[keep, "target_multiherb_fraction"].to_numpy(float),
    ])
    X = np.column_stack([np.ones(X.shape[0]), X])
    yy = y[keep]
    xtx_inv = np.linalg.pinv(X.T @ X)
    beta = xtx_inv @ X.T @ yy
    residual = yy - X @ beta
    leverage = np.clip(np.diag(X @ xtx_inv @ X.T), 0, 1 - 1e-12)
    scaled = residual / (1 - leverage)
    meat = X.T @ ((scaled ** 2)[:, None] * X)
    cov_hc3 = xtx_inv @ meat @ xtx_inv
    se_hc3 = np.sqrt(np.clip(np.diag(cov_hc3), 0, None))
    lofo_rows.append({
        "omitted_formula_id": row.formula_id,
        "omitted_formula_name_en": row.formula_name_en,
        "n_remaining": int(keep.sum()),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "hc3_multiherb_coefficient": float(beta[-1]),
        "hc3_multiherb_se": float(se_hc3[-1]),
    })
lofo = pd.DataFrame(lofo_rows)
lofo.to_csv(ROB / "cross_formula_leave_one_out_sensitivity.csv", index=False, encoding="utf-8-sig")

rng = np.random.default_rng(20260802)
boot = []
for rep in range(10000):
    idx = rng.integers(0, len(benchmark), size=len(benchmark))
    if len(np.unique(x[idx])) < 2 or len(np.unique(y[idx])) < 2:
        continue
    rho, _ = spearmanr(x[idx], y[idx])
    boot.append(rho)
boot = np.asarray(boot, float)
pd.DataFrame([{
    "n_formulae": len(benchmark),
    "spearman_rho": float(rho_full),
    "spearman_p": float(p_full),
    "bootstrap_replicates_requested": 10000,
    "bootstrap_replicates_valid": len(boot),
    "bootstrap_q025": float(np.quantile(boot, 0.025)),
    "bootstrap_q975": float(np.quantile(boot, 0.975)),
    "lofo_rho_min": float(lofo.spearman_rho.min()),
    "lofo_rho_max": float(lofo.spearman_rho.max()),
    "seed": 20260802,
    "interpretive_boundary": "Metric-behavior check under union coverage; not an independent biological association.",
}]).to_csv(ROB / "cross_formula_association_sensitivity_summary.csv", index=False, encoding="utf-8-sig")


# Explicit SQRCD rank, denominator and percentile within the prespecified comparison band.
sq = benchmark.loc[benchmark.formula_id.eq("F000")].iloc[0]
band = benchmark.loc[benchmark.comparison_band.eq(sq.comparison_band)].sort_values("target_fdi", ascending=True).reset_index(drop=True)
rank_low_to_high = int(band.index[band.formula_id.eq("F000")][0]) + 1
rank_high_to_low = len(band) - rank_low_to_high + 1
pd.DataFrame([{
    "formula_id": "F000",
    "comparison_band": sq.comparison_band,
    "n_in_band": len(band),
    "rank_low_to_high_fdi": rank_low_to_high,
    "rank_high_to_low_fdi": rank_high_to_low,
    "percentile_within_band": float(sq.target_fdi_percentile_within_band),
    "target_fdi": float(sq.target_fdi),
}]).to_csv(ROB / "sqrcd_within_band_rank.csv", index=False, encoding="utf-8-sig")


def make_figure6() -> Path:
    df = benchmark.sort_values("target_fdi", ascending=True).reset_index(drop=True)
    palette = {"compound": "#D97706", "target": "#2563EB", "sqrcd": "#C026D3", "ink": "#1F2937", "grid": "#D1D5DB"}
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 10.2), constrained_layout=True)
    ax = axes[0, 0]
    ax.axis("off")
    boxes = [(0.04, "20\nprespecified", "#E0E7FF"), (0.39, "14\nbenchmark retained", "#DBEAFE"), (0.74, "15\nwith SQRCD", "#EDE9FE")]
    for x0, txt, color in boxes:
        ax.add_patch(plt.Rectangle((x0, .34), .22, .34, transform=ax.transAxes, facecolor=color, edgecolor=palette["ink"], lw=1.2))
        ax.text(x0+.11, .51, txt, transform=ax.transAxes, ha="center", va="center", fontsize=12, weight="bold", color=palette["ink"])
    for x1, x2 in [(.27,.39),(.62,.74)]:
        ax.annotate("", xy=(x2,.51), xytext=(x1,.51), xycoords=ax.transAxes, arrowprops=dict(arrowstyle="->", color=palette["ink"], lw=1.6))
    ax.text(.5,.13,"Identity gate: government formulary-HERB 2.0 Jaccard >=0.80; herb mapping >=0.80",transform=ax.transAxes,ha="center",fontsize=9,color="#4B5563")
    ax.set_title("A  Prespecified portability stress test", loc="left", fontsize=13, weight="bold")

    ax = axes[0, 1]
    yy = np.arange(len(df))
    ax.hlines(yy, df.target_fdi, df.compound_fdi, color="#CBD5E1", lw=2)
    ax.scatter(df.compound_fdi, yy, color=palette["compound"], s=42, label="Compound FDI", zorder=3)
    ax.scatter(df.target_fdi, yy, color=palette["target"], s=42, label="Target FDI", zorder=3)
    ax.axvline(.5, color="#111827", ls="--", lw=1)
    ax.text(.498, len(df)-1.35, "non-redundant upper bound", ha="right", va="top", fontsize=8, color="#374151")
    ax.set_yticks(yy, df.formula_name_en, fontsize=8)
    ax.set_xlim(0, .515)
    ax.set_xlabel("Formula Discrimination Index (0-0.5)")
    ax.grid(axis="x", color=palette["grid"], lw=.6)
    ax.legend(frameon=False, ncol=2, loc="lower left", fontsize=8)
    ax.set_title("B  Compound and target FDI, sorted by target FDI", loc="left", fontsize=13, weight="bold")

    ax = axes[1, 0]
    ax.scatter(benchmark.target_multiherb_fraction, benchmark.target_fdi, s=45, color=palette["target"], edgecolor="white", lw=.7)
    labels = set(["F000", benchmark.loc[benchmark.target_fdi.idxmin(),"formula_id"], benchmark.loc[benchmark.target_fdi.idxmax(),"formula_id"]])
    for _, r in benchmark.loc[benchmark.formula_id.isin(labels)].iterrows():
        ax.annotate(r.formula_id, (r.target_multiherb_fraction, r.target_fdi), xytext=(5,5), textcoords="offset points", fontsize=8, weight="bold" if r.formula_id=="F000" else "normal")
    ax.set_ylim(0, .5)
    ax.set_xlabel("Cross-herb target redundancy fraction")
    ax.set_ylabel("Target FDI (0-0.5)")
    ax.grid(color=palette["grid"], lw=.6)
    ax.text(.02,.03,f"Spearman rho={rho_full:.3f}; LOFO range {lofo.spearman_rho.min():.3f} to {lofo.spearman_rho.max():.3f}",transform=ax.transAxes,fontsize=8,color="#4B5563")
    ax.set_title("C  Expected metric behavior under union coverage", loc="left", fontsize=13, weight="bold")

    ax = axes[1, 1]
    bands = [x for x in ["4-6", "7-10", "11-20"] if x in set(benchmark.comparison_band)]
    for j, band_name in enumerate(bands):
        temp = benchmark.loc[benchmark.comparison_band.eq(band_name)]
        offsets = np.linspace(-.12,.12,len(temp)) if len(temp)>1 else np.array([0.0])
        colors = [palette["sqrcd"] if f=="F000" else palette["target"] for f in temp.formula_id]
        ax.scatter(np.full(len(temp),j)+offsets, temp.target_fdi, s=[75 if f=="F000" else 42 for f in temp.formula_id], c=colors, edgecolor="white", lw=.7)
        if "F000" in set(temp.formula_id):
            rr=temp.loc[temp.formula_id.eq("F000")].iloc[0]
            ax.annotate("SQRCD",(j,rr.target_fdi),xytext=(-8,6),ha="right",textcoords="offset points",fontsize=8,weight="bold",color=palette["sqrcd"])
    ax.axhline(.5,color="#111827",ls="--",lw=1)
    ax.set_xticks(range(len(bands)),[f"{b}\n(n={int((benchmark.comparison_band==b).sum())})" for b in bands])
    ax.set_ylim(0,.515)
    ax.set_ylabel("Target FDI (0-0.5)")
    ax.grid(axis="y",color=palette["grid"],lw=.6)
    ax.set_title("D  Target FDI by formula-size comparison band", loc="left", fontsize=13, weight="bold")
    fig.suptitle("Cross-formula portability benchmark\nStress testing EQF across independently specified formulae", fontsize=15, weight="bold", color=palette["ink"])
    out = EXT / "figures" / "Figure6_external_benchmark_V7.png"
    fig.savefig(out, dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def make_compound_sensitivity_figure() -> Path:
    plot = entity.loc[entity.status.eq("estimated")].pivot(index=["formula_id","formula_name_en"], columns="entity_representation", values="compound_fdi").reset_index()
    plot = plot.sort_values("normalized_compound_name")
    yv = np.arange(len(plot))
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.2, 7.2), gridspec_kw={"width_ratios":[1.65,1]}, constrained_layout=True)
    ax.hlines(yv, plot.normalized_compound_name, plot.ETCM_native_record_id, color="#94A3B8", lw=2)
    ax.scatter(plot.ETCM_native_record_id, yv, label="ETCM native record ID", color="#D97706", s=42, zorder=3)
    ax.scatter(plot.normalized_compound_name, yv, label="Normalized compound name", color="#2563EB", s=42, zorder=3)
    ax.axvline(.5, color="#111827", ls="--", lw=1)
    ax.set_yticks(yv, plot.formula_name_en, fontsize=8)
    ax.set_xlim(.475,.501)
    ax.set_xlabel("Compound FDI (focused 0.475-0.500 scale)")
    ax.grid(axis="x", color="#D1D5DB", lw=.6)
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    ax.set_title("A  FDI under two available entity rules", loc="left", fontsize=12, weight="bold")
    delta=(plot.ETCM_native_record_id-plot.normalized_compound_name)*1000
    ax2.barh(yv,delta,color="#2563EB",edgecolor="#1E3A8A",height=.55)
    ax2.axvline(0,color="#111827",lw=.8)
    ax2.set_yticks(yv,[])
    ax2.set_xlabel("FDI decrease after name normalization (x10^-3)")
    ax2.grid(axis="x",color="#D1D5DB",lw=.6)
    ax2.set_title("B  Change attributable to name normalization", loc="left", fontsize=12, weight="bold")
    fig.suptitle("Compound entity harmonization sensitivity", fontsize=15, weight="bold")
    fig.text(.5,.005,"PubChem CID and parent InChIKey representations were not estimable because uniform identifiers were absent.",ha="center",fontsize=8,color="#4B5563")
    out=FIG/"Supplementary_Figure_S3_compound_entity_harmonization.png"
    fig.savefig(out,dpi=320,bbox_inches="tight",facecolor="white")
    fig.savefig(out.with_suffix(".pdf"),bbox_inches="tight",facecolor="white")
    plt.close(fig)
    return out


fig6 = make_figure6()
figs3 = make_compound_sensitivity_figure()

summary_json = {
    "rater_A": rater_a_path.name,
    "rater_B": rater_b_path.name,
    "blank_template": blank_path.name,
    "four_state_raw_agreement": float((merged.rater_state_A == merged.rater_state_B).mean()),
    "four_state_kappa": float(agreement_summary.loc[(agreement_summary.analysis.eq("four_state")) & (agreement_summary.metric.eq("cohen_kappa_unweighted")),"value"].iloc[0]),
    "collapsed_raw_agreement": float((collapsed_a == collapsed_b).mean()),
    "collapsed_kappa": float(agreement_summary.loc[(agreement_summary.analysis.eq("collapsed_unavailable")) & (agreement_summary.metric.eq("cohen_kappa_unweighted")),"value"].iloc[0]),
    "disagreements": len(disagreements),
    "stopped_vs_not_estimable_disagreements": int(((merged.rater_state_A.eq("stopped") & merged.rater_state_B.eq("not estimable")) | (merged.rater_state_B.eq("stopped") & merged.rater_state_A.eq("not estimable"))).sum()),
    "compound_native_fdi_range": [float(entity.loc[(entity.entity_representation.eq("ETCM_native_record_id")) & entity.status.eq("estimated"),"compound_fdi"].min()), float(entity.loc[(entity.entity_representation.eq("ETCM_native_record_id")) & entity.status.eq("estimated"),"compound_fdi"].max())],
    "compound_normalized_name_fdi_range": [float(entity.loc[(entity.entity_representation.eq("normalized_compound_name")) & entity.status.eq("estimated"),"compound_fdi"].min()), float(entity.loc[(entity.entity_representation.eq("normalized_compound_name")) & entity.status.eq("estimated"),"compound_fdi"].max())],
    "pubchem_parent_inchikey_status": "not estimable",
    "rho_full": float(rho_full),
    "rho_lofo_range": [float(lofo.spearman_rho.min()), float(lofo.spearman_rho.max())],
    "sqrcd_band_n": len(band),
    "sqrcd_rank_low_to_high": rank_low_to_high,
    "sqrcd_percentile": float(sq.target_fdi_percentile_within_band),
    "figure6": str(fig6),
    "supplementary_figure_s3": str(figs3),
    "independence_boundary": "Distinct initials and dates are present, but rater independence cannot be established from the workbooks alone.",
}
(ROOT / "logs" / "V6_analysis_summary.json").write_text(json.dumps(summary_json, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(summary_json, ensure_ascii=False, indent=2))
