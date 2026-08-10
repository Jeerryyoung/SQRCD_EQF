from __future__ import annotations

import hashlib
import json
import re
import time
from pathlib import Path

import pandas as pd
import requests


BASE = Path(__file__).resolve().parents[1]
API = "http://www.tcmip.cn:18124"
HERB_FORMULA_URL = "http://47.92.70.12/static/download_data/V2/HERB_formula_info_v2.txt"
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "SQRCD-EQF-V5/1.0 (academic reproducibility audit)"})


def norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def split_herbs(value: str) -> list[str]:
    return [x.strip() for x in str(value).split("|") if x.strip()]


def split_herb2(value: str) -> list[str]:
    return [re.sub(r"\s+", "", x.strip()) for x in str(value).split(",") if x.strip() and x.strip() != "NA"]


def get_json(method: str, path: str, **kwargs) -> dict:
    last = None
    for attempt in range(4):
        try:
            response = SESSION.request(method, API + path, timeout=60, **kwargs)
            response.raise_for_status()
            payload = response.json()
            if payload.get("code") != 1:
                raise RuntimeError(payload.get("msg", "ETCM response code was not 1"))
            return payload
        except Exception as exc:  # network retry is logged by caller output files
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"ETCM request failed after retries: {path}: {last}")


def exact_formula_records(query: str) -> tuple[list[dict], dict]:
    search = get_json(
        "POST",
        "/home/homeEsSearch/",
        json={
            "search_text": query.lower(),
            "type": "traditional_chinese_medicine_formula",
            "pageNo": 1,
            "pageSize": 100,
        },
    )
    uuids = search.get("result", {}).get("name_list", [])
    records: list[dict] = []
    raw = {"search": search, "browse": []}
    for uid in uuids:
        page = get_json(
            "POST",
            "/home/browse/",
            json={
                "type": "traditional_chinese_medicine_formula",
                "name_list": uid,
                "pageNo": 1,
                "pageSize": 1000,
            },
        )
        raw["browse"].append(page)
        for block in page.get("data", []):
            for row in block.get("data", []):
                names = row.get("Formula Name in Chinese Pinyin", [])
                name = names[0] if isinstance(names, list) and names else str(names)
                if norm(name) == norm(query):
                    records.append(row)
    return records, raw


def alias_map() -> dict[str, str]:
    aliases = pd.read_csv(BASE / "config" / "herb_aliases.csv")
    return {norm(row.input_name): row.etcm_name for row in aliases.itertuples()}


def harmonize(herbs: list[str], aliases: dict[str, str]) -> list[str]:
    return sorted({aliases.get(norm(h), h) for h in herbs}, key=norm)


def jaccard(left: list[str], right: list[str]) -> float:
    a, b = {norm(x) for x in left}, {norm(x) for x in right}
    return len(a & b) / len(a | b) if (a or b) else 0.0


def herb_profile(herb: str) -> dict:
    return get_json("GET", "/home/detail/", params={"id": herb, "type": "herb"})


def parse_profile(herb: str, payload: dict) -> tuple[list[dict], list[dict], dict]:
    compounds: dict[str, dict] = {}
    targets: dict[str, dict] = {}
    metadata = {"herb": herb, "latin_name": "", "medicinal_part": ""}
    for section in payload.get("data", []):
        if section.get("key") == "Basic Information":
            for item in section.get("value", []):
                if item.get("key") == "Herb Name in Latin":
                    metadata["latin_name"] = re.sub(r"<[^>]+>", "", str(item.get("value", "")))
                if item.get("key") == "Medicinal Part":
                    metadata["medicinal_part"] = item.get("value", "")
        if section.get("key") == "Network Visualization":
            for node in section.get("value", {}).get("nodes", []):
                node_id = str(node.get("id", ""))
                if node_id.startswith("TCMIP-I-"):
                    compounds[node_id] = {
                        "canonical_herb_id": herb,
                        "node_id": node_id,
                        "node_label": node.get("label", ""),
                        "source": "ETCM_v2_herb_network",
                    }
        if section.get("key") == "Component Target":
            for relation in section.get("value", {}).get("data", []):
                genes = relation.get("Gene Symbol", [])
                ingredients = relation.get("Ingredient Name", [])
                genes = genes if isinstance(genes, list) else [genes]
                ingredients = ingredients if isinstance(ingredients, list) else [ingredients]
                for gene in genes:
                    if not gene:
                        continue
                    key = str(gene).upper()
                    targets[key] = {
                        "canonical_herb_id": herb,
                        "node_id": key,
                        "node_label": key,
                        "source": "ETCM_v2_component_target",
                    }
                # Recover ingredient names not represented in the network block.
                for ingredient in ingredients:
                    if ingredient and not any(x["node_label"] == ingredient for x in compounds.values()):
                        pseudo = "NAME-" + hashlib.sha1(str(ingredient).encode("utf-8")).hexdigest()[:12]
                        compounds[pseudo] = {
                            "canonical_herb_id": herb,
                            "node_id": pseudo,
                            "node_label": ingredient,
                            "source": "ETCM_v2_component_target_name",
                        }
    return list(compounds.values()), list(targets.values()), metadata


def main() -> None:
    candidates = pd.read_csv(BASE / "config" / "formula_candidates.csv")
    aliases = alias_map()
    raw_formula_dir = BASE / "raw" / "formulas"
    raw_herb_dir = BASE / "raw" / "herb_profiles"
    processed = BASE / "processed"
    for directory in [raw_formula_dir, raw_herb_dir, processed]:
        directory.mkdir(parents=True, exist_ok=True)

    herb_formula_path = BASE / "raw" / "HERB_formula_info_v2.txt"
    if not herb_formula_path.exists():
        response = SESSION.get(HERB_FORMULA_URL, timeout=120)
        response.raise_for_status()
        herb_formula_path.write_bytes(response.content)
    herb2 = pd.read_csv(herb_formula_path, sep="\t", dtype=str).fillna("")
    herb2["_name_norm"] = herb2["Formula_pinyin_name"].map(norm)

    audit_rows, locked_rows, exclusion_rows = [], [], []
    for row in candidates.itertuples():
        primary = harmonize(split_herbs(row.moh_herbs), aliases)
        if row.formula_id == "F000":
            selected, records, raw = primary, [], {"investigator_locked": True}
            concordance, etcm_concordance = 1.0, 1.0
            source_record = "investigator_locked_SQRCD"
        else:
            records, raw = exact_formula_records(row.etcm_query)
            etcm_scored = []
            for record in records:
                herbs = record.get("Herbs Contained in This Formula (Chinese Pinyin)", [])
                herbs = harmonize(herbs if isinstance(herbs, list) else [herbs], aliases)
                etcm_scored.append((jaccard(primary, herbs), herbs, record))
            if etcm_scored:
                etcm_concordance, etcm_selected, etcm_best = max(etcm_scored, key=lambda x: (x[0], -abs(len(x[1]) - len(primary))))
            else:
                etcm_concordance, etcm_selected, etcm_best = 0.0, [], {}
            herb2_rows = herb2[herb2["_name_norm"] == norm(row.etcm_query)]
            herb2_scored = []
            for _, hrow in herb2_rows.iterrows():
                herbs = harmonize(split_herb2(hrow.get("Herbs_in_pinyin", "")), aliases)
                herb2_scored.append((jaccard(primary, herbs), herbs, hrow.to_dict()))
            if herb2_scored:
                concordance, selected, best = max(herb2_scored, key=lambda x: (x[0], -abs(len(x[1]) - len(primary))))
                source_record = f"HERB2:{best.get('Formula_id','')}:{best.get('Source','')}"
            else:
                concordance, selected, source_record = 0.0, [], "no_exact_HERB2_record"
            raw["HERB2_exact_records"] = [x[2] for x in herb2_scored]
            raw["ETCM_best_herbs"] = etcm_selected
        (raw_formula_dir / f"{row.formula_id}_{row.etcm_query}.json").write_text(
            json.dumps(raw, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        included = len(primary) >= 4 and len(primary) <= 25 and concordance >= 0.80
        audit_rows.append(
            {
                "formula_id": row.formula_id,
                "formula_name_en": row.formula_name_en,
                "formula_name_zh": row.formula_name_zh,
                "primary_n": len(primary),
                "etcm_n": len(selected),
                "cross_source_jaccard": concordance,
                "etcm_best_jaccard": etcm_concordance,
                "exact_etcm_records": len(records),
                "exact_herb2_records": 0 if row.formula_id == "F000" else len(herb2_scored),
                "selected_etcm_source": source_record,
                "included_identity_gate": included,
                "primary_herbs": "|".join(primary),
                "selected_secondary_herbs": "|".join(selected),
            }
        )
        if not included:
            exclusion_rows.append(
                {
                    "formula_id": row.formula_id,
                    "formula_name_en": row.formula_name_en,
                    "stage": "identity_concordance",
                    "reason": "cross-source Jaccard below 0.80 or formula size outside 4-25",
                    "observed_jaccard": concordance,
                }
            )
            continue
        for herb in primary:
            locked_rows.append(
                {
                    "formula_id": row.formula_id,
                    "formula_name_en": row.formula_name_en,
                    "formula_name_zh": row.formula_name_zh,
                    "indication_group": row.indication_group,
                    "canonical_herb_id": herb,
                    "identity_source_1": row.source_primary,
                    "identity_source_2": "ETCM_v2" if row.formula_id == "F000" else "HERB_2.0",
                    "cross_source_jaccard": concordance,
                }
            )

    audit = pd.DataFrame(audit_rows)
    locked = pd.DataFrame(locked_rows)
    exclusions = pd.DataFrame(exclusion_rows)
    audit.to_csv(processed / "formula_identity_audit.csv", index=False, encoding="utf-8-sig")

    compound_edges, target_edges, herb_meta, failed_herbs = [], [], [], []
    for herb in sorted(locked["canonical_herb_id"].unique(), key=norm):
        try:
            cache_path = raw_herb_dir / f"{herb}.json"
            if cache_path.exists():
                payload = json.loads(cache_path.read_text(encoding="utf-8"))
            else:
                payload = herb_profile(herb)
                cache_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
            compounds, targets, metadata = parse_profile(herb, payload)
            compound_edges.extend(compounds)
            target_edges.extend(targets)
            metadata.update({"compound_count": len(compounds), "target_count": len(targets), "mapped": bool(compounds or targets)})
            herb_meta.append(metadata)
        except Exception as exc:
            failed_herbs.append({"canonical_herb_id": herb, "error": str(exc)})
            herb_meta.append({"herb": herb, "latin_name": "", "medicinal_part": "", "compound_count": 0, "target_count": 0, "mapped": False})

    meta = pd.DataFrame(herb_meta)
    mapped = set(meta.loc[meta["mapped"], "herb"])
    if failed_herbs:
        pd.DataFrame(failed_herbs).to_csv(processed / "herb_retrieval_failures.csv", index=False, encoding="utf-8-sig")
    mapping_rows, retained_formula_ids = [], []
    for formula_id, group in locked.groupby("formula_id"):
        herbs = sorted(group["canonical_herb_id"].unique(), key=norm)
        mapped_n = sum(h in mapped for h in herbs)
        rate = mapped_n / len(herbs)
        keep = rate >= 0.80
        mapping_rows.append({"formula_id": formula_id, "herb_n": len(herbs), "mapped_n": mapped_n, "mapping_rate": rate, "included_mapping_gate": keep})
        if keep:
            retained_formula_ids.append(formula_id)
        else:
            exclusions.loc[len(exclusions)] = [formula_id, group.iloc[0]["formula_name_en"], "herb_mapping", "ETCM herb mapping rate below 0.80", rate]

    locked = locked[locked["formula_id"].isin(retained_formula_ids)].copy()
    locked.to_csv(processed / "formulas_locked.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(mapping_rows).to_csv(processed / "herb_mapping_audit.csv", index=False, encoding="utf-8-sig")
    exclusions.to_csv(processed / "benchmark_exclusions.csv", index=False, encoding="utf-8-sig")
    meta.to_csv(processed / "herb_metadata.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(compound_edges).drop_duplicates(["canonical_herb_id", "node_id"]).to_csv(processed / "herb_compound_edges.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(target_edges).drop_duplicates(["canonical_herb_id", "node_id"]).to_csv(processed / "herb_target_edges.csv", index=False, encoding="utf-8-sig")
    print(json.dumps({"candidate_formulae": len(candidates), "identity_pass": int(audit["included_identity_gate"].sum()), "mapping_pass": len(retained_formula_ids), "unique_herbs": int(locked["canonical_herb_id"].nunique()), "compound_edges": len(compound_edges), "target_edges": len(target_edges)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
