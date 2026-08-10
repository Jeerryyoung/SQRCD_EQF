from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "05_figures"
EXT_FIG = ROOT / "external_benchmark" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
EXT_FIG.mkdir(parents=True, exist_ok=True)

INK = "#17324D"
BLUE = "#2F6B8A"
BLUE_LIGHT = "#DCEAF2"
GOLD = "#C7922F"
GOLD_LIGHT = "#F4E8C8"
ORANGE = "#D46A32"
OLIVE = "#657A3D"
GREY = "#68737D"
GREY_LIGHT = "#E9EDF0"
PINK = "#B85C7A"


def save(fig, path: Path) -> None:
    fig.savefig(path.with_suffix(".png"), dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def figure1_framework() -> None:
    fig, ax = plt.subplots(figsize=(15.5, 7.8))
    ax.set_xlim(0, 15.5)
    ax.set_ylim(0, 8)
    ax.axis("off")

    steps = [
        (0.25, "Evidence\ndossier", "Identity, chemistry, targets,\nhuman context and provenance", BLUE),
        (3.25, "Gate-specific\nassignment", "Qualified | partially qualified |\nnot estimable | stopped", OLIVE),
        (6.25, "Independent\nassessment", "Two blinded field-level records;\nconclusions withheld", GOLD),
        (9.25, "Reproducibility\naudit", "Agreement, kappa, domain map\nand disagreement taxonomy", ORANGE),
        (12.25, "Claim-eligibility\nstate", "Concordant state retained;\nboundaries remain explicit", PINK),
    ]
    for x, title, body, color in steps:
        box = FancyBboxPatch(
            (x, 4.38), 2.45, 1.78,
            boxstyle="round,pad=0.05,rounding_size=0.12",
            linewidth=1.5, edgecolor=color, facecolor="white"
        )
        ax.add_patch(box)
        ax.text(x + 1.225, 5.72, title, ha="center", va="center", fontsize=9.5, weight="bold", color=color, linespacing=1.0)
        ax.text(x + 1.225, 4.92, body, ha="center", va="center", fontsize=7.9, color=INK, linespacing=1.2)
    for x in [2.70, 5.70, 8.70, 11.70]:
        ax.add_patch(FancyArrowPatch((x, 5.27), (x + 0.5, 5.27), arrowstyle="-|>", mutation_scale=14, linewidth=1.4, color=GREY))

    ax.text(0.25, 7.4, "Evidence Qualification Framework (EQF)", fontsize=18, weight="bold", color=INK)
    ax.text(0.25, 7.0, "Evidence-state assignment, reproducibility and cross-formula portability before mechanistic interpretation", fontsize=10.8, color=GREY)

    gate_y = 2.45
    gates = [
        ("Formula\nidentity", "mandatory"),
        ("Preparation\nchemistry", "mandatory"),
        ("Chemical\nreconstruction", "mandatory"),
        ("Target\nprovenance", "mandatory"),
        ("Human\ncontext", "supportive"),
        ("Genetics", "supportive"),
        ("Perturbation", "supportive"),
        ("Structural\neligibility", "dependent"),
    ]
    colors = {"mandatory": BLUE, "supportive": GOLD, "dependent": PINK}
    for i, (label, group) in enumerate(gates):
        x = 0.35 + i * 1.86
        ax.add_patch(FancyBboxPatch((x, gate_y), 1.56, 0.92, boxstyle="round,pad=0.03,rounding_size=0.08", facecolor=colors[group], edgecolor="none", alpha=0.92))
        ax.text(x + 0.78, gate_y + 0.46, label, ha="center", va="center", color="white", fontsize=8.2, weight="bold")
    ax.text(0.35, 3.62, "Eight prespecified evidence domains", fontsize=11.2, weight="bold", color=INK)

    ax.text(0.35, 1.58, "Framework output", fontsize=10.5, weight="bold", color=INK)
    ax.text(2.65, 1.58, "Evidence-state profile", fontsize=9.3, color=BLUE)
    ax.text(5.15, 1.58, "Decision trace", fontsize=9.3, color=OLIVE)
    ax.text(7.45, 1.58, "Reproducibility record", fontsize=9.3, color=ORANGE)
    ax.text(10.65, 1.58, "Explicit stopping point", fontsize=9.3, color=PINK)
    ax.text(0.35, 0.68, "Independent assessment evaluates reproducibility of state assignment, not biological correctness.", fontsize=10, style="italic", color=GREY)
    save(fig, FIG / "Figure1_EQF_evidence_state_framework_V7")


def lofo_figure() -> None:
    df = pd.read_csv(ROOT / "02_robustness" / "cross_formula_leave_one_out_sensitivity.csv").sort_values("spearman_rho")
    full = -0.9857142857142858
    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    y = np.arange(len(df))
    ax.scatter(df["spearman_rho"], y, s=46, color=BLUE, edgecolor="white", linewidth=0.7, zorder=3)
    ax.axvline(full, color=ORANGE, linewidth=1.7, linestyle="--", label=f"All formulae: rho={full:.3f}")
    ax.set_yticks(y)
    ax.set_yticklabels(df["omitted_formula_name_en"], fontsize=8.2)
    ax.set_xlim(-0.990, -0.978)
    ax.set_xlabel("Spearman rho after omitting one formula")
    fig.subplots_adjust(top=0.84)
    fig.suptitle("Leave-one-formula-out metric-behavior sensitivity", x=0.16, y=0.98, ha="left", fontsize=15, weight="bold", color=INK)
    fig.text(0.16, 0.925, "15 omissions; the redundancy-FDI relation remains definition-linked", fontsize=9.3, color=GREY)
    ax.grid(axis="x", color=GREY_LIGHT, linewidth=0.8)
    ax.legend(frameon=False, loc="lower right", fontsize=8.8)
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, FIG / "Supplementary_Figure_S4_LOFO_portability_V7")


def domain_portability_heatmap() -> None:
    df = pd.read_csv(ROOT / "external_benchmark" / "results" / "benchmark_gate_profiles.csv")
    columns = ["identity", "chemical_topology", "target_topology", "preparation_chemistry", "human_context", "genetics", "perturbation", "structural_eligibility"]
    labels = ["Identity", "Chemical\ntopology", "Target\ntopology", "Preparation\nchemistry", "Human\ncontext", "Genetics", "Perturbation", "Structural\neligibility"]
    state_code = {"qualified": 0, "partially qualified": 1, "not estimable": 2, "stopped": 3}
    arr = df[columns].replace(state_code).to_numpy(dtype=float)
    cmap = ListedColormap([BLUE, GOLD, GREY_LIGHT, PINK])
    fig, ax = plt.subplots(figsize=(10.5, 7.1))
    ax.imshow(arr, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5)
    ax.set_xticks(range(len(columns)), labels, fontsize=8.4)
    ax.set_yticks(range(len(df)), df["formula_name_en"], fontsize=8.2)
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, pad=5)
    for i in range(arr.shape[0] + 1):
        ax.axhline(i - 0.5, color="white", linewidth=0.8)
    for j in range(arr.shape[1] + 1):
        ax.axvline(j - 0.5, color="white", linewidth=0.8)
    fig.subplots_adjust(top=0.78, bottom=0.13)
    fig.suptitle("EQF domain-state portability across independently specified formulae", x=0.16, y=0.98, ha="left", fontsize=15, weight="bold", color=INK)
    fig.text(0.16, 0.925, "Topology modules were portable; absent downstream evidence remained explicitly unavailable", fontsize=9.2, color=GREY)
    legend = [Patch(facecolor=c, label=l) for c, l in zip([BLUE, GOLD, GREY_LIGHT, PINK], ["Qualified", "Partially qualified", "Not estimable", "Stopped"])]
    ax.legend(handles=legend, ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.06), loc="upper center", fontsize=8.6)
    ax.spines[:].set_visible(False)
    save(fig, FIG / "Supplementary_Figure_S5_domain_portability_V7")


def inter_rater_map() -> None:
    df = pd.read_csv(ROOT / "02_eqf" / "second_rater_field_level_joined.csv")
    formula_order = df[["dossier_id", "formula_name_en_A"]].drop_duplicates().sort_values("dossier_id")
    item_order = df[["rubric_code", "domain_A"]].drop_duplicates()
    domain_order = ["Identity", "Preparation chemistry", "Chemical reconstruction", "Target provenance", "Human context", "Genetics", "Perturbation", "Structural eligibility"]
    item_order["domain_order"] = item_order["domain_A"].map({d: i for i, d in enumerate(domain_order)})
    item_order["prefix"] = item_order["rubric_code"].str.extract(r"(\d+)").astype(int)
    item_order = item_order.sort_values(["domain_order", "prefix"])
    f_index = {x: i for i, x in enumerate(formula_order["dossier_id"])}
    r_index = {x: i for i, x in enumerate(item_order["rubric_code"])}
    mat = np.zeros((len(formula_order), len(item_order)), dtype=int)
    for row in df.itertuples():
        i = f_index[row.dossier_id]
        j = r_index[row.rubric_code]
        if row.rater_state_A == row.rater_state_B:
            value = 0
        elif {row.rater_state_A, row.rater_state_B} == {"stopped", "not estimable"}:
            value = 1
        else:
            value = 2
        mat[i, j] = value

    fig = plt.figure(figsize=(14.5, 7.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[5.2, 1.25], wspace=0.3)
    ax = fig.add_subplot(gs[0, 0])
    cmap = ListedColormap([BLUE_LIGHT, ORANGE, PINK])
    ax.imshow(mat, aspect="auto", cmap=cmap, vmin=-0.5, vmax=2.5)
    ax.set_yticks(range(len(formula_order)), formula_order["formula_name_en_A"], fontsize=8.1)
    ax.set_xticks(range(len(item_order)), item_order["rubric_code"], rotation=90, fontsize=6.6)
    starts = []
    for domain in domain_order:
        idx = np.where(item_order["domain_A"].to_numpy() == domain)[0]
        if len(idx):
            starts.append((idx[0], idx[-1], domain))
    for start, end, domain in starts:
        ax.axvline(start - 0.5, color="white", linewidth=1.4)
        ax.text((start + end) / 2, -1.0, domain.replace(" ", "\n"), ha="center", va="bottom", fontsize=7.1, color=INK)
    ax.axvline(mat.shape[1] - 0.5, color="white", linewidth=1.4)
    ax.set_xlabel("50 prespecified rubric fields")
    ax.spines[:].set_visible(False)

    ax2 = fig.add_subplot(gs[0, 1])
    labels = ["Exact agreement", "Stopped vs\nnot estimable", "Other"]
    values = [689, 60, 1]
    colors = [BLUE, ORANGE, PINK]
    bars = ax2.barh(range(3), values, color=colors, edgecolor="white")
    ax2.set_yticks(range(3), labels, fontsize=8.5)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 750)
    ax2.set_xlabel("Ratings")
    ax2.set_title("Disagreement\ntaxonomy", loc="left", fontsize=11.7, weight="bold", color=INK)
    for bar, value in zip(bars, values):
        ax2.text(value + 8, bar.get_y() + bar.get_height()/2, str(value), va="center", fontsize=9, color=INK)
    ax2.grid(axis="x", color=GREY_LIGHT, linewidth=0.8)
    ax2.spines[["top", "right"]].set_visible(False)
    fig.subplots_adjust(top=0.76, bottom=0.12)
    fig.suptitle("Evidence-state reproducibility across 15 formulae and 50 rubric fields", x=0.05, y=0.98, ha="left", fontsize=15.5, weight="bold", color=INK)
    fig.text(0.05, 0.925, "Four-state agreement: 689/750; 60 of 61 disagreements localize to one downstream-state boundary", fontsize=9.5, color=GREY)
    save(fig, FIG / "Supplementary_Figure_S6_inter_rater_disagreement_V7")

    taxonomy = pd.DataFrame({
        "disagreement_type": ["stopped versus not estimable", "partially qualified versus qualified"],
        "n": [60, 1],
        "share_of_disagreements": [60 / 61, 1 / 61],
        "interpretation": [
            "Systematic downstream-state definition boundary across structural fields",
            "Single SQRCD PE04 evidence-threshold disagreement",
        ],
    })
    taxonomy.to_csv(ROOT / "02_eqf" / "second_rater_disagreement_taxonomy.csv", index=False)


if __name__ == "__main__":
    figure1_framework()
    lofo_figure()
    domain_portability_heatmap()
    inter_rater_map()
    print("V7 figures and disagreement taxonomy created")
