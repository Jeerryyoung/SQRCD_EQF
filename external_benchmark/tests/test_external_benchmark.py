from pathlib import Path
import gzip

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parents[1]
P = BASE / "processed"
R = BASE / "results"


def main() -> None:
    formulas = pd.read_csv(P / "formulas_locked.csv")
    identity = pd.read_csv(P / "formula_identity_audit.csv")
    mapping = pd.read_csv(P / "herb_mapping_audit.csv")
    summary = pd.read_csv(R / "external_benchmark_summary.csv")
    shapley = pd.read_csv(R / "formula_shapley.csv")
    gates = pd.read_csv(R / "benchmark_gate_profiles.csv")

    assert summary["formula_id"].nunique() == 15
    assert (summary["formula_id"] != "F000").sum() == 14
    assert formulas.groupby("formula_id")["canonical_herb_id"].nunique().between(4, 25).all()
    assert mapping.loc[mapping["formula_id"].isin(summary["formula_id"]), "mapping_rate"].ge(0.80).all()
    assert identity.loc[identity["formula_id"].isin(summary["formula_id"]), "cross_source_jaccard"].ge(0.80).all()
    assert summary["compound_fdi"].between(0, 0.5 + 1e-12).all()
    assert summary["target_fdi"].between(0, 0.5 + 1e-12).all()
    assert (summary["target_fdi"] < summary["compound_fdi"]).all()
    efficiency = shapley.groupby(["formula_id", "node_type"])["shapley"].sum()
    assert np.allclose(efficiency, 1.0)
    expected_rows = 2 * sum((2 ** n) - 1 for n in formulas.groupby("formula_id")["canonical_herb_id"].nunique())
    with gzip.open(R / "complete_subformula_lattice.csv.gz", "rt", encoding="utf-8") as handle:
        observed_rows = sum(1 for _ in handle) - 1
    assert observed_rows == expected_rows == 257410
    assert (gates["structural_eligibility"] == "stopped").all()
    assert not (gates.drop(columns=["formula_id", "formula_name_en"]) == "second-rater complete").any().any()
    print("External benchmark invariants passed across identity, mapping, topology, attribution, enumeration and stopping rules.")


if __name__ == "__main__":
    main()
