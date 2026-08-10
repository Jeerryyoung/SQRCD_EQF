from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    summary = pd.read_csv(ROOT / "external_benchmark/results/external_benchmark_summary.csv")
    gates = pd.read_csv(ROOT / "external_benchmark/results/benchmark_gate_profiles.csv")
    agreement = pd.read_csv(ROOT / "results/eqf/second_rater_agreement_summary.csv")
    taxonomy = pd.read_csv(ROOT / "results/eqf/second_rater_disagreement_taxonomy.csv")
    crate = json.loads((ROOT / "ro-crate-metadata.json").read_text(encoding="utf-8"))

    assert summary["formula_id"].nunique() == 15
    assert (summary["formula_id"] != "F000").sum() == 14
    assert summary["compound_fdi"].between(0, 0.5 + 1e-12).all()
    assert summary["target_fdi"].between(0, 0.5 + 1e-12).all()
    assert (summary["target_fdi"] < summary["compound_fdi"]).all()
    assert (gates["structural_eligibility"] == "stopped").all()

    values = agreement.set_index(["analysis", "metric"])["value"]
    assert np.isclose(values.loc[("four_state", "raw_agreement")], 689 / 750)
    assert int(taxonomy["n"].sum()) == 61
    assert crate["@graph"][1]["codeRepository"].endswith("/SQRCD_EQF")
    assert not any(ROOT.rglob("*.xlsx"))
    print("Public-release invariants passed.")


if __name__ == "__main__":
    main()
