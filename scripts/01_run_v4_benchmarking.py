"""Run V4 formula-specificity benchmarking analyses.

Outputs quantify annotation topology. They do not estimate efficacy, synergy,
clinical indispensability, exposure, or dose response.
"""
from __future__ import annotations

from collections import Counter
from itertools import combinations
from pathlib import Path
import hashlib, json, math

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

V4=Path(__file__).resolve().parents[1]; SCHEME=V4.parent; PROJECT=SCHEME.parent
REBUILD=PROJECT/"SQRCD_rebuild_2026-07-28"; RNG=np.random.default_rng(20260801)
HERBS=["Nanshashen","Baizhu","Fuling","Gancao","Chenpi","Banxia","Chaihu","Baishao","Zhike","Huangqi","Beishashen","Doukou","Yiyiren","Kuxingren","Yuliren","Baiziren"]
DOSE=dict(zip(HERBS,[30,30,15,6,15,15,15,30,15,60,30,10,20,10,30,30]))
TOTAL_MASS=sum(DOSE.values())

def truthy(s): return s.astype(str).str.strip().str.lower().isin({"true","1","yes","y"})
def union(sets): return set().union(*sets) if sets else set()
def gini(x):
    x=np.sort(np.asarray(x,dtype=float)); n=len(x)
    if n==0 or x.sum()==0:return np.nan
    return float((2*np.sum((np.arange(n)+1)*x)/(n*x.sum()))-(n+1)/n)
def entropy_metrics(x):
    x=np.asarray(x,dtype=float); x=x[x>0]; p=x/x.sum(); h=float(-(p*np.log(p)).sum())
    return h/math.log(len(x)),math.exp(h),float(np.sum(p*p)),gini(x)
def curve_from_support(counts,n=16):
    """Exact mean retention over all subsets, including the empty endpoint."""
    counts=np.asarray(counts,dtype=int); out=[]
    for removed in range(n+1):
        retained=[]
        for r in counts:
            miss=math.comb(n-r,removed)/math.comb(n,removed) if removed<=n-r else 0.0
            retained.append(1-miss)
        out.append(np.mean(retained) if len(retained) else 0.0)
    return np.asarray(out)
def fdi_from_support(counts,n=16):
    curve=curve_from_support(counts,n); x=np.arange(n+1)/n
    return float(1-np.trapezoid(curve,x)),curve
def bitset(items,index):
    z=0
    for x in items:z|=1<<index[x]
    return z
def shapley_union(sets,universe):
    support=Counter(x for s in sets.values() for x in s); n=max(len(universe),1)
    return {h:sum(1/support[x] for x in sets[h])/n for h in HERBS},support
def concentration(vals):
    x=np.asarray(list(vals),dtype=float); p=x/x.sum()
    return float(np.sum(p*p)),gini(p),float(np.max(p))

def main():
    for d in ["01_lineage","02_eqf","02_robustness","03_context","05_figures","06_supplementary","07_manuscript","08_reproducibility","logs"]:(V4/d).mkdir(exist_ok=True)
    chem_path=REBUILD/"03_constituent_evidence"/"constituent_evidence_all_source_rows.csv"
    rel_path=SCHEME/"09_tables"/"experimental_compound_target_disease_relations.csv"
    dis_path=REBUILD/"04_disease_target_evidence"/"disease_gene_evidence_matrix.csv"
    identity_path=V4/"00_protocol"/"formula_identity_lock.csv"
    chem=pd.read_csv(chem_path,low_memory=False); rel=pd.read_csv(rel_path,low_memory=False); dis=pd.read_csv(dis_path,low_memory=False)
    chem=chem[chem.herb_pinyin.isin(HERBS)].copy(); rel=rel[rel.herb_pinyin.isin(HERBS)].copy()
    chem["compound_id"]=chem.evidence_compound_key.fillna(chem.compound_key).astype(str)
    rel["gene"]=rel.gene.astype(str).str.strip(); rel["disease_supported_bool"]=truthy(rel.disease_supported)
    dr=rel[rel.disease_supported_bool & rel.gene.ne("")].copy()
    hc={h:set(chem.loc[chem.herb_pinyin.eq(h),"compound_id"]) for h in HERBS}
    hg={h:set(dr.loc[dr.herb_pinyin.eq(h),"gene"]) for h in HERBS}
    fc=union(list(hc.values())); fg=union(list(hg.values()))

    # V4 gate profile: categorical primary output, no additive score in main inference.
    gates=pd.DataFrame([
      [1,"Formula identity","mandatory","Exact membership, species, processing, dose","qualified","Mechanistic workflow may proceed"],
      [2,"Preparation chemistry","mandatory","Preparation-matched measured constituents","partially qualified","Exploratory chemistry only"],
      [3,"Target provenance","mandatory","Traceable assay class and human target","partially qualified","Annotation coverage only"],
      [4,"Human disease context","supportive","Phenotype- and tissue-matched human data","partially qualified","Disease context only"],
      [5,"Human genetics","supportive","Reproducible variant-level and regional inputs","not estimable","No causal promotion"],
      [6,"Formula perturbation","supportive","Formula-specific intervention signature","not estimable","No reversal ranking"],
      [7,"Machine-learning discrimination","supportive","Independent labelled formulae and grouped/external validation","not estimable","Supervised task not fitted"],
      [8,"Structural eligibility","dependent","Qualified compound-protein pair surviving mandatory gates","stopped","No docking or MD"],
    ],columns=["layer","evidence_domain","gate_role","required_input","state","claim_eligibility"])
    gates.to_csv(V4/"02_eqf"/"V4_EQF_gate_profile.csv",index=False,encoding="utf-8-sig")
    pd.DataFrame([{"state":s,"n_domains":int((gates.state==s).sum())} for s in ["qualified","partially qualified","not estimable","stopped"]]).to_csv(V4/"02_eqf"/"V4_gate_state_counts.csv",index=False,encoding="utf-8-sig")

    # Exact enumeration of every non-empty subformula.
    gi={g:i for i,g in enumerate(sorted(fg))}; ci={c:i for i,c in enumerate(sorted(fc))}
    gb=[bitset(hg[h],gi) for h in HERBS]; cb=[bitset(hc[h],ci) for h in HERBS]
    gene_bits=[0]*65536; comp_bits=[0]*65536; rows=[]
    gene_edge_total=sum(DOSE[h]*len(hg[h]) for h in HERBS); comp_edge_total=sum(DOSE[h]*len(hc[h]) for h in HERBS)
    for mask in range(1,65536):
        lsb=mask & -mask; i=lsb.bit_length()-1; prev=mask^lsb
        gene_bits[mask]=gene_bits[prev]|gb[i]; comp_bits[mask]=comp_bits[prev]|cb[i]
        inds=[j for j in range(16) if mask>>j&1]; n=len(inds)
        g=gene_bits[mask].bit_count(); c=comp_bits[mask].bit_count()
        mass=sum(DOSE[HERBS[j]] for j in inds)
        gw=sum(DOSE[HERBS[j]]*len(hg[HERBS[j]]) for j in inds)/gene_edge_total
        cw=sum(DOSE[HERBS[j]]*len(hc[HERBS[j]]) for j in inds)/comp_edge_total
        rows.append([mask,n,16-n,"|".join(HERBS[j] for j in inds),mass,mass/TOTAL_MASS,g,g/len(fg),c,c/len(fc),gw,cw])
    subsets=pd.DataFrame(rows,columns=["subset_mask","n_retained","n_removed","retained_herbs","crude_mass_g","crude_mass_fraction","disease_genes","target_retention","unique_compounds","compound_retention","crude_mass_weighted_target_edge_load","crude_mass_weighted_compound_edge_load"])
    subsets.to_csv(V4/"02_robustness"/"all_65535_nonempty_subformulae.csv",index=False,encoding="utf-8-sig")
    qrows=[]
    for n,x in subsets.groupby("n_retained"):
        qrows.append({"n_retained":n,"n_removed":16-n,"n_combinations":len(x),
          **{f"target_{k}":v for k,v in zip(["mean","median","q05","q95","min","max"],[x.target_retention.mean(),x.target_retention.median(),x.target_retention.quantile(.05),x.target_retention.quantile(.95),x.target_retention.min(),x.target_retention.max()])},
          **{f"compound_{k}":v for k,v in zip(["mean","median","q05","q95","min","max"],[x.compound_retention.mean(),x.compound_retention.median(),x.compound_retention.quantile(.05),x.compound_retention.quantile(.95),x.compound_retention.min(),x.compound_retention.max()])},
          "fraction_target_retention_ge_0_80":float((x.target_retention>=.8).mean()),
          "mass_weighted_target_edge_load_mean":x.crude_mass_weighted_target_edge_load.mean(),
          "mass_weighted_compound_edge_load_mean":x.crude_mass_weighted_compound_edge_load.mean()})
    curves=pd.DataFrame(qrows).sort_values("n_retained"); curves.to_csv(V4/"02_robustness"/"complete_subformula_retention_summary.csv",index=False,encoding="utf-8-sig")

    # Transparent Formula Discrimination Index.
    gene_shap,g_support=shapley_union(hg,fg); comp_shap,c_support=shapley_union(hc,fc)
    target_fdi,target_curve=fdi_from_support(list(g_support.values())); compound_fdi,compound_curve=fdi_from_support(list(c_support.values()))
    fdi=pd.DataFrame([
      ["disease-supported targets",target_fdi,1-target_fdi,len(fg),"unique-node retention"],
      ["provisional compounds",compound_fdi,1-compound_fdi,len(fc),"unique-node retention"],
    ],columns=["annotation_unit","FDI","normalized_retention_AUC","full_nodes","interpretation"])
    fdi.to_csv(V4/"02_robustness"/"formula_discrimination_index.csv",index=False,encoding="utf-8-sig")
    pd.DataFrame({"n_removed":range(17),"fraction_removed":np.arange(17)/16,"mean_target_retention":target_curve,"mean_compound_retention":compound_curve}).to_csv(V4/"02_robustness"/"exact_mean_retention_curve.csv",index=False,encoding="utf-8-sig")

    # Exact Shapley values for union coverage; unique and redundant components.
    shap=[]
    for h in HERBS:
        ug=sum(g_support[g]==1 for g in hg[h])/len(fg); uc=sum(c_support[c]==1 for c in hc[h])/len(fc)
        shap.append([h,DOSE[h],gene_shap[h],ug,gene_shap[h]-ug,comp_shap[h],uc,comp_shap[h]-uc,len(hg[h]),len(hc[h])])
    shap=pd.DataFrame(shap,columns=["herb","crude_mass_g","target_shapley","target_unique_component","target_redundant_component","compound_shapley","compound_unique_component","compound_redundant_component","herb_target_degree","herb_compound_degree"])
    shap.to_csv(V4/"02_robustness"/"exact_shapley_annotation_contributions.csv",index=False,encoding="utf-8-sig")

    # Direct redundancy and concentration metrics.
    ge,geff,ghhi,ggini=entropy_metrics(list(g_support.values())); ce,ceff,chhi,cgini=entropy_metrics(list(c_support.values()))
    target_sc=concentration(gene_shap.values()); comp_sc=concentration(comp_shap.values())
    metrics=pd.DataFrame([
      ["target",len(fg),ge,geff,ghhi,ggini,np.mean(np.array(list(g_support.values()))>=2),np.mean(list(g_support.values())),*target_sc],
      ["compound",len(fc),ce,ceff,chhi,cgini,np.mean(np.array(list(c_support.values()))>=2),np.mean(list(c_support.values())),*comp_sc],
    ],columns=["unit","nodes","normalized_entropy","entropy_equivalent_equally_supported_nodes","support_HHI","support_Gini","multiherb_support_fraction","mean_herb_support","shapley_HHI","shapley_Gini","maximum_shapley_share"])
    metrics.to_csv(V4/"02_robustness"/"redundancy_concentration_metrics.csv",index=False,encoding="utf-8-sig")

    # Evidence-tier target composition and tier-specific FDI.
    tiers={
      "all_disease_supported":dr,
      "ChEMBL_direct_quantitative_le_10uM":dr[dr.relation_source.eq("ChEMBL_human_single_protein_le_10uM")],
      "TCMSP_labelled_validated":dr[dr.relation_source.eq("TCMSP_labelled_validated")],
      "human_genetic_disease_support":dr[truthy(dr.human_genetic_evidence)],
      "curated_or_clinical_disease_support":dr[truthy(dr.curated_or_clinical_support)],
      "GeneCards_score_ge_3_support":dr[truthy(dr.genecards_score_ge_3)],
    }
    tier_rows=[]
    for name,df in tiers.items():
        sets={h:set(df.loc[df.herb_pinyin.eq(h),"gene"]) for h in HERBS}; universe=union(list(sets.values()))
        if universe:
            _,sup=shapley_union(sets,universe); val,_=fdi_from_support(list(sup.values())); multi=float(np.mean(np.array(list(sup.values()))>=2))
        else: val=np.nan; multi=np.nan
        tier_rows.append([name,len(df),df[["herb_pinyin","gene"]].drop_duplicates().shape[0],df.MOL_ID.nunique(),len(universe),val,multi])
    pd.DataFrame(tier_rows,columns=["evidence_tier","relation_rows","unique_herb_gene_edges","unique_compounds","unique_genes","target_FDI","multiherb_target_fraction"]).to_csv(V4/"02_robustness"/"target_evidence_tier_sensitivity.csv",index=False,encoding="utf-8-sig")

    # Disease evidence formation chain for the 139 genes.
    dis=dis[dis.gene_symbol.astype(str).isin(fg)].copy()
    source_defs={
      "Open Targets primary phenotype score >0":pd.to_numeric(dis.opentargets_primary_score,errors="coerce").fillna(0)>0,
      "Open Targets genetic evidence score >0":pd.to_numeric(dis.opentargets_genetic_score,errors="coerce").fillna(0)>0,
      "Open Targets clinical evidence score >0":pd.to_numeric(dis.opentargets_clinical_score,errors="coerce").fillna(0)>0,
      "GWAS Catalog association present":pd.to_numeric(dis.gwas_catalog_associations,errors="coerce").fillna(0)>0,
      "GeneCards relevance score >=3":pd.to_numeric(dis.genecards_relevance_score,errors="coerce").fillna(0)>=3,
      "OMIM gene-map row present":pd.to_numeric(dis.omim_gene_map_rows,errors="coerce").fillna(0)>0,
    }
    pd.DataFrame([[k,int(v.sum()),"supportive annotation; not automatically causal"] for k,v in source_defs.items()],columns=["disease_evidence_source_rule","unique_genes","interpretation"]).to_csv(V4/"06_supplementary"/"disease_target_evidence_source_composition.csv",index=False,encoding="utf-8-sig")
    dis.to_csv(V4/"06_supplementary"/"disease_supported_139_gene_evidence_matrix.csv",index=False,encoding="utf-8-sig")

    # Annotation perturbation: independent edge thinning rather than resampling with replacement.
    orig_loss=np.array([len(fg-union([hg[x] for x in HERBS if x!=h])) for h in HERBS]); orig_top=set(np.array(HERBS)[np.argsort(-orig_loss)[:3]])
    pert=[]
    for rate in [.9,.8,.7,.6]:
        for rep in range(1000):
            sets={h:set(x for x in hg[h] if RNG.random()<rate) for h in HERBS}; uf=union(list(sets.values()))
            loss=np.array([len(uf-union([sets[x] for x in HERBS if x!=h])) for h in HERBS]); top=set(np.array(HERBS)[np.argsort(-loss)[:3]])
            sup=Counter(x for s in sets.values() for x in s); f,_=fdi_from_support(list(sup.values())) if sup else (np.nan,None)
            tau=kendalltau(orig_loss,loss,nan_policy="omit").statistic
            pert.append([rate,rep+1,len(uf),f,tau,len(top&orig_top)/len(top|orig_top),int(np.argmax(loss)!=np.argmax(orig_loss))])
    pert=pd.DataFrame(pert,columns=["edge_retention_rate","replicate","target_coverage","target_FDI","herb_loss_rank_tau","top3_jaccard","top1_rank_reversal"])
    pert.to_csv(V4/"02_robustness"/"annotation_edge_thinning_replicates.csv",index=False,encoding="utf-8-sig")
    summ=[]
    for rate,x in pert.groupby("edge_retention_rate"):
        row={"edge_retention_rate":rate,"replicates":len(x),"top1_rank_reversal_probability":x.top1_rank_reversal.mean()}
        for col in ["target_coverage","target_FDI","herb_loss_rank_tau","top3_jaccard"]:
            row.update({f"{col}_median":x[col].median(),f"{col}_q025":x[col].quantile(.025),f"{col}_q975":x[col].quantile(.975)})
        summ.append(row)
    pd.DataFrame(summ).to_csv(V4/"02_robustness"/"annotation_edge_thinning_summary.csv",index=False,encoding="utf-8-sig")

    # Degree-preserving bipartite rewiring: appropriate for attribution, not FDI.
    edges={(i,g) for i,h in enumerate(HERBS) for g in hg[h]}; edge_list=list(edges); null=[]
    obs_hhi=target_sc[0]; obs_max=max(gene_shap.values()); obs_loss_hhi=np.sum((orig_loss/orig_loss.sum())**2) if orig_loss.sum() else np.nan
    for rep in range(1000):
        for _ in range(len(edge_list)*3):
            a,b=RNG.integers(0,len(edge_list),2); h1,g1=edge_list[a]; h2,g2=edge_list[b]
            if h1==h2 or g1==g2 or (h1,g2) in edges or (h2,g1) in edges: continue
            edges.remove((h1,g1)); edges.remove((h2,g2)); edges.add((h1,g2)); edges.add((h2,g1)); edge_list[a]=(h1,g2); edge_list[b]=(h2,g1)
        sets={HERBS[i]:set(g for h,g in edges if h==i) for i in range(16)}; uf=union(list(sets.values()))
        sv,_=shapley_union(sets,uf); hhi,_,mx=concentration(sv.values())
        losses=np.array([len(uf-union([sets[x] for x in HERBS if x!=h])) for h in HERBS]); lhhi=np.sum((losses/losses.sum())**2) if losses.sum() else np.nan
        null.append([rep+1,hhi,mx,lhhi])
    null=pd.DataFrame(null,columns=["replicate","shapley_HHI","maximum_shapley_share","LOHO_loss_HHI"]); null.to_csv(V4/"02_robustness"/"degree_preserving_rewiring_null.csv",index=False,encoding="utf-8-sig")
    pd.DataFrame([
      ["shapley_HHI",obs_hhi,null.shapley_HHI.mean(),null.shapley_HHI.std(),(1+(null.shapley_HHI>=obs_hhi).sum())/(len(null)+1)],
      ["maximum_shapley_share",obs_max,null.maximum_shapley_share.mean(),null.maximum_shapley_share.std(),(1+(null.maximum_shapley_share>=obs_max).sum())/(len(null)+1)],
      ["LOHO_loss_HHI",obs_loss_hhi,null.LOHO_loss_HHI.mean(),null.LOHO_loss_HHI.std(),(1+(null.LOHO_loss_HHI>=obs_loss_hhi).sum())/(len(null)+1)],
    ],columns=["metric","observed","null_mean","null_sd","empirical_upper_tail_p"]).to_csv(V4/"02_robustness"/"degree_preserving_null_summary.csv",index=False,encoding="utf-8-sig")

    # Pareto frontier over herb count, target retention and compound retention.
    pareto=[]
    for n,x in subsets.groupby("n_retained"):
        y=x.sort_values(["target_retention","compound_retention"],ascending=False)
        best=-1
        for _,r in y.iterrows():
            if r.compound_retention>best:
                pareto.append(r); best=r.compound_retention
    pd.DataFrame(pareto).to_csv(V4/"02_robustness"/"subformula_pareto_frontier.csv",index=False,encoding="utf-8-sig")

    # External validation and inter-rater agreement are audited, not fabricated.
    pd.DataFrame([
      ["Independent formula identities","No harmonized external formula benchmark package supplied","available in principle","identity-only benchmarking insufficient"],
      ["External herb-constituent universe","ITCMDB chemistry available, but cross-database standardized formula-specific extraction not executed","partial","cannot support target-level comparison"],
      ["External herb-target evidence","No uniformly processed target-evidence universe beyond the 16 SQRCD herbs","not estimable","do not claim cross-formula validation"],
      ["External disease context","Indication-specific disease evidence not harmonized across formulae","not estimable","do not rank external formulae"],
      ["Final transferability claim","SQRCD case-tested framework","case tested","proposed framework; not externally validated"],
    ],columns=["requirement","observed","state","decision"]).to_csv(V4/"06_supplementary"/"external_formula_benchmark_feasibility_audit.csv",index=False,encoding="utf-8-sig")
    pd.DataFrame([[d,"",""] for d in gates.evidence_domain],columns=["evidence_domain","reviewer_2_state","adjudication_note"]).to_csv(V4/"06_supplementary"/"independent_rater_template.csv",index=False,encoding="utf-8-sig")

    # Data lineage and deterministic summary.
    paths=[chem_path,rel_path,dis_path,identity_path]
    pd.DataFrame([[str(p),p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()] for p in paths],columns=["path","bytes","sha256"]).to_csv(V4/"01_lineage"/"V4_source_manifest.csv",index=False,encoding="utf-8-sig")
    summary={"version":"V4","seed":20260801,"herbs":16,"daily_crude_mass_g":TOTAL_MASS,"unique_compounds":len(fc),"disease_supported_genes":len(fg),"nonempty_subformulae":len(subsets),"target_FDI":target_fdi,"compound_FDI":compound_fdi,"edge_thinning_replicates":len(pert),"degree_preserving_null_replicates":len(null),"external_formula_validation":"not_estimable"}
    (V4/"logs"/"V4_run_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2))

if __name__=="__main__": main()
