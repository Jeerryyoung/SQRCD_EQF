from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import ListedColormap

V4=Path(__file__).resolve().parents[1]; F=V4/"05_figures"; F.mkdir(exist_ok=True)
NAVY="#18324B"; BLUE="#2F6B9A"; TEAL="#2A9D8F"; GOLD="#E9A23B"; RED="#C8553D"; PLUM="#7A5195"; GREY="#747C89"; LIGHT="#EEF3F6"; DARK="#263238"
plt.rcParams.update({"font.family":"Arial","font.size":9,"axes.titlesize":12,"axes.labelsize":9.5,"axes.spines.top":False,"axes.spines.right":False})
def save(fig,name): fig.savefig(F/name,dpi=320,bbox_inches="tight",facecolor="white"); plt.close(fig)
def panel(ax,letter,title,sub=""):
    if letter: ax.text(-.06,1.08,letter,transform=ax.transAxes,fontweight="bold",fontsize=13)
    ax.set_title(title,loc="left",fontweight="bold",color=NAVY,pad=12)
    if sub: ax.text(0,1.01,sub,transform=ax.transAxes,color=GREY,fontsize=8)

# Figure 1: framework and algorithmic gate roles.
g=pd.read_csv(V4/"02_eqf"/"V4_EQF_gate_profile.csv")
fig,ax=plt.subplots(figsize=(11,6)); ax.axis("off")
ax.text(.03,.94,"Evidence Qualification Framework",fontsize=18,fontweight="bold",color=NAVY)
ax.text(.03,.89,"Gate-specific claim eligibility takes precedence over an additive score",fontsize=10,color=GREY)
cols={"qualified":TEAL,"partially qualified":GOLD,"not estimable":GREY,"stopped":RED}
for i,r in g.iterrows():
    x=.04+(i%4)*.24; y=.60-(i//4)*.27; col=cols[r.state]
    ax.add_patch(FancyBboxPatch((x,y),.205,.19,boxstyle="round,pad=.012",fc="white",ec=col,lw=2))
    ax.text(x+.102,y+.145,r.evidence_domain,ha="center",va="center",fontweight="bold",fontsize=9)
    ax.text(x+.102,y+.08,r.gate_role.upper(),ha="center",fontsize=7,color=NAVY)
    ax.text(x+.102,y+.025,r.state.upper(),ha="center",fontsize=7,color="white",bbox=dict(boxstyle="round,pad=.22",fc=col,ec="none"))
ax.add_patch(FancyBboxPatch((.04,.08),.92,.12,boxstyle="round,pad=.012",fc=LIGHT,ec="none"))
ax.text(.07,.15,"Algorithm",fontweight="bold",color=NAVY); ax.text(.17,.15,"lock identity → qualify chemistry and provenance → test context and uncertainty → assign claim eligibility → stop or proceed",fontsize=9)
ax.text(.07,.105,"Case output",fontweight="bold",color=NAVY); ax.text(.20,.105,"1 qualified · 3 partial · 3 not estimable · 1 stopped",fontsize=9)
save(fig,"Figure1_EQF_gate_algorithm.png")

# Figure 2: source and evidence tier sensitivity.
src=pd.read_csv(V4/"02_robustness"/"source_database_ablation_summary.csv")
tier=pd.read_csv(V4/"02_robustness"/"target_evidence_tier_sensitivity.csv")
ds=pd.read_csv(V4/"06_supplementary"/"disease_target_evidence_source_composition.csv")
fig,axs=plt.subplots(1,3,figsize=(13,4.8),gridspec_kw={"width_ratios":[1,1.15,1]})
panel(axs[0],"a","Chemical source sensitivity","Provisional unique compounds")
axs[0].barh(range(len(src)),src.unique_compounds,color=[BLUE if n==1 else TEAL for n in src.n_databases]); axs[0].set_yticks(range(len(src)),src.database_set); axs[0].invert_yaxis(); axs[0].set_xlabel("Compounds")
panel(axs[1],"b","Target-evidence tier sensitivity","FDI and eligible genes")
y=np.arange(len(tier)); axs[1].barh(y,tier.target_FDI,color=PLUM); axs[1].set_yticks(y,tier.evidence_tier.str.replace("_"," ")); axs[1].invert_yaxis(); axs[1].set_xlim(0,.52); axs[1].set_xlabel("Target FDI")
for i,r in tier.iterrows(): axs[1].text(r.target_FDI+.01,i,f"n={int(r.unique_genes)}",va="center",fontsize=7)
panel(axs[2],"c","Disease-evidence composition","Genes among the 139-gene set")
axs[2].barh(range(len(ds)),ds.unique_genes,color=GOLD); axs[2].set_yticks(range(len(ds)),ds.disease_evidence_source_rule); axs[2].invert_yaxis(); axs[2].set_xlabel("Unique genes")
for i,v in enumerate(ds.unique_genes): axs[2].text(v+1,i,str(v),va="center",fontsize=8)
fig.tight_layout(w_pad=3); save(fig,"Figure2_source_evidence_sensitivity.png")

# Figure 3: context compatibility and STC QC.
qc=pd.read_csv(V4/"03_context"/"GSE245885_QC_summary.csv"); loo=pd.read_csv(V4/"03_context"/"GSE245885_leave_one_sample_out_QC.csv")
ctx=pd.DataFrame([
 ["GSE245885","STC","colon","subject","No","Cancer-surgery\ncontext","Low-confidence\ncontext"],
 ["GSE36701","IBS-C","rectal mucosa","repeated biopsy","No","Participant\nblocking","Disease context"],
 ["GSE166869","IBS-C","duodenum/\njejunum","subject","No","Tissue-specific\nnull","Sensitivity only"],
],columns=["Dataset","Phenotype","Tissue","Unit","Intervention","Critical feature","EQF use"])
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,4.7),gridspec_kw={"width_ratios":[1.45,1]})
panel(ax1,"a","Phenotype–tissue–intervention compatibility","Disease expression is not a formula perturbation")
ax1.axis("off"); tbl=ax1.table(cellText=ctx.values,colLabels=ctx.columns,cellLoc="center",loc="center",bbox=[0,0,1,.88]); tbl.auto_set_font_size(False); tbl.set_fontsize(6.8)
for (r,c),cell in tbl.get_celld().items():
    cell.set_edgecolor("white"); cell.set_facecolor(NAVY if r==0 else ("#F7E5C1" if r==1 else LIGHT)); cell.get_text().set_color("white" if r==0 else DARK); cell.get_text().set_fontweight("bold" if r==0 else "normal")
panel(ax2,"b","GSE245885 influence audit","Broad differential expression persists after single-sample deletion")
ax2.bar(np.arange(len(loo)),loo.FDR05_abs_logFC1,color=[BLUE if g=="NP" else RED for g in loo.omitted_group]); ax2.axhline(qc.full_FDR05_abs_logFC1.iloc[0],color=DARK,ls="--",label="Full model")
ax2.set_xticks(range(len(loo)),loo.omitted_sample,rotation=45,ha="right"); ax2.set_ylabel("DE genes: FDR<0.05 and |log2FC|≥1"); ax2.legend(frameon=False)
ax2.text(.02,.98,f"MDS label-permutation P={qc.MDS_centroid_permutation_p.iloc[0]:.3f}\nDesign remains clinically confounded",transform=ax2.transAxes,va="top",fontsize=8,color=RED)
fig.tight_layout(w_pad=3); save(fig,"Figure3_context_compatibility.png")

# Figure 4: all subsets and FDI.
cur=pd.read_csv(V4/"02_robustness"/"complete_subformula_retention_summary.csv").sort_values("n_removed")
sub=pd.read_csv(V4/"02_robustness"/"all_65535_nonempty_subformulae.csv",usecols=["n_retained","crude_mass_fraction","target_retention"])
fdi=pd.read_csv(V4/"02_robustness"/"formula_discrimination_index.csv")
fig,axs=plt.subplots(1,3,figsize=(13,4.6))
panel(axs[0],"a","Complete subformula retention","All 65,535 non-empty subsets")
x=cur.n_removed.to_numpy(); axs[0].fill_between(x,cur.target_q05,cur.target_q95,color=BLUE,alpha=.18); axs[0].plot(x,cur.target_median,"o-",color=BLUE,label="Targets")
axs[0].fill_between(x,cur.compound_q05,cur.compound_q95,color=GOLD,alpha=.18); axs[0].plot(x,cur.compound_median,"s-",color=GOLD,label="Compounds"); axs[0].set_xlabel("Herbs removed"); axs[0].set_ylabel("Retention from F00"); axs[0].set_ylim(0,1.03); axs[0].legend(frameon=False)
panel(axs[1],"b","Crude-mass-weighted sensitivity","Crude mass is not exposure")
hb=axs[1].hexbin(sub.crude_mass_fraction,sub.target_retention,gridsize=35,cmap="Blues",mincnt=1,bins="log"); axs[1].set_xlabel("Retained crude-mass fraction"); axs[1].set_ylabel("Target retention"); fig.colorbar(hb,ax=axs[1],label="log10 subset count")
panel(axs[2],"c","Formula Discrimination Index","FDI = 1 − normalized retention AUC")
axs[2].bar([0,1],fdi.FDI,color=[BLUE,GOLD]); axs[2].set_xticks([0,1],["Targets","Compounds"]); axs[2].set_ylim(0,.55); axs[2].set_ylabel("FDI (higher = faster annotation loss)")
for i,v in enumerate(fdi.FDI): axs[2].text(i,v+.015,f"{v:.3f}",ha="center",fontweight="bold")
fig.tight_layout(w_pad=3); save(fig,"Figure4_complete_dropout_FDI.png")

# Figure 5: perturbation, Shapley and null.
pt=pd.read_csv(V4/"02_robustness"/"annotation_edge_thinning_summary.csv").sort_values("edge_retention_rate")
sh=pd.read_csv(V4/"02_robustness"/"exact_shapley_annotation_contributions.csv").sort_values("target_shapley")
dn=pd.read_csv(V4/"02_robustness"/"degree_preserving_null_summary.csv")
fig,axs=plt.subplots(1,3,figsize=(13.2,5.2),gridspec_kw={"width_ratios":[1,1.25,1]})
panel(axs[0],"a","Annotation perturbation stability","1,000 replicates per edge-retention rate")
x=pt.edge_retention_rate; med=pt.herb_loss_rank_tau_median; lo=pt.herb_loss_rank_tau_q025; hi=pt.herb_loss_rank_tau_q975
axs[0].errorbar(x,med,yerr=[med-lo,hi-med],fmt="o-",color=BLUE,capsize=4); axs[0].set_xlim(.55,.95); axs[0].set_ylim(0,1.02); axs[0].set_xlabel("Retained herb–target edges"); axs[0].set_ylabel("Herb-loss rank Kendall τ")
ax=axs[0].twinx(); ax.plot(x,pt.top1_rank_reversal_probability,"s--",color=RED); ax.set_ylabel("Top-rank reversal probability",color=RED); ax.set_ylim(0,.55)
panel(axs[1],"b","Exact Shapley annotation contribution","Union-coverage game; values sum to 1")
y=np.arange(len(sh)); axs[1].barh(y,sh.target_unique_component,color=TEAL,label="Unique"); axs[1].barh(y,sh.target_redundant_component,left=sh.target_unique_component,color=GOLD,label="Redundant"); axs[1].set_yticks(y,sh.herb); axs[1].set_xlabel("Target Shapley contribution"); axs[1].legend(frameon=False,fontsize=8)
panel(axs[2],"c","Degree-preserving rewiring null","Attribution concentration, not FDI")
yy=np.arange(len(dn)); axs[2].scatter(dn.observed,yy,color=RED,label="Observed",zorder=3); axs[2].errorbar(dn.null_mean,yy,xerr=1.96*dn.null_sd,fmt="o",color=GREY,label="Null mean ±1.96 SD"); axs[2].set_yticks(yy,dn.metric.str.replace("_"," ")); axs[2].set_xlabel("Concentration metric"); axs[2].legend(frameon=False,fontsize=8)
fig.tight_layout(w_pad=3); save(fig,"Figure5_perturbation_shapley_null.png")

# Figure 6: transferability audit and claim map.
ext=pd.read_csv(V4/"06_supplementary"/"external_formula_benchmark_feasibility_audit.csv")
fig,ax=plt.subplots(figsize=(11,5.8)); ax.axis("off")
ax.text(.03,.94,"Claim-eligibility and transferability audit",fontsize=17,fontweight="bold",color=NAVY)
ax.text(.03,.89,"The framework is case-tested; cross-formula validation requires a harmonized external herb–target universe",fontsize=10,color=GREY)
items=[("Qualified","Identity-locked composition and deterministic rebuilding",TEAL),("Exploratory","Source/tier sensitivity, FDI, Shapley and annotation perturbation",GOLD),("Not estimable","Causal genetics, formula perturbation and supervised discrimination",GREY),("Stopped","Structural simulation without an eligible compound–protein pair",RED)]
for i,(a,b,c) in enumerate(items):
    y=.69-i*.14; ax.add_patch(FancyBboxPatch((.04,y),.18,.09,boxstyle="round,pad=.01",fc=c,ec="none")); ax.text(.13,y+.045,a,ha="center",va="center",color="white",fontweight="bold")
    ax.add_patch(FancyBboxPatch((.25,y),.70,.09,boxstyle="round,pad=.01",fc=LIGHT,ec="none")); ax.text(.28,y+.045,b,va="center",fontsize=9)
ax.text(.04,.08,"Transferability status",fontweight="bold",color=NAVY); ax.text(.20,.08,"Proposed and case-tested; external benchmark not estimable from the supplied target universe.",fontsize=10)
save(fig,"Figure6_claim_transferability_map.png")
print("V4 figures built")
