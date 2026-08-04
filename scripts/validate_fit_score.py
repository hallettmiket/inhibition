import sys, re, glob, os
sys.path.insert(0,'/home/UWO/twu383/repos/inhibition')
import numpy as np, pandas as pd
from shared import pose_vector as pv

RB="/data/lab_vm/append_only/inhibition/05_redock_benchmark"
REC="/data/lab_vm/immutable/inhibition/receptor/6VAJ_prepared.pdb"

def pdb_heavy(path, want_resi=None):
    """residue -> (n,3) heavy atoms, from a PDB."""
    out={}
    for ln in open(path):
        if not ln.startswith(("ATOM","HETATM")): continue
        el=ln[76:78].strip() or ln[12:16].strip()[0]
        if el.upper()=="H": continue
        try: r=int(ln[22:26]); xyz=(float(ln[30:38]),float(ln[38:46]),float(ln[46:54]))
        except ValueError: continue
        if want_resi is not None and r not in want_resi: continue
        out.setdefault(r,[]).append(xyz)
    return {k:np.array(v) for k,v in out.items()}

def lig_heavy_pdb(path):
    pts=[]
    for ln in open(path):
        if not ln.startswith(("ATOM","HETATM")): continue
        el=(ln[76:78].strip() or ln[12:16].strip()[0]).upper()
        if el=="H": continue
        pts.append((float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
    return np.array(pts)

def pdbqt_models(path):
    models,cur=[],[]
    for ln in open(path):
        if ln.startswith("MODEL"): cur=[]
        elif ln.startswith("ENDMDL"):
            if cur: models.append(np.array(cur))
            cur=[]
        elif ln.startswith(("ATOM","HETATM")):
            el=ln[77:79].strip().upper() if len(ln)>78 else ""
            if el.startswith("H") and el!="HG": continue
            try: cur.append((float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
            except ValueError: pass
    if cur: models.append(np.array(cur))
    return models

sys.path.insert(0,'/home/UWO/twu383/repos/inhibition/integration/app')
import pose3d as p3d
RESI=tuple(p3d.pocket_resi())
rec=pdb_heavy(REC, set(RESI))
print(f"receptor pocket residues resolved: {len(rec)}/{len(RESI)}")

cases=[]
for ref in sorted(glob.glob(f"{RB}/cases_1/refs_6vaj/*_ref6vaj.pdb")):
    cid=os.path.basename(ref).replace("_ref6vaj.pdb","")
    dock=f"{RB}/dock_1/cross/poses/{cid}_out.pdbqt"
    if not os.path.exists(dock): continue
    try:
        cx=lig_heavy_pdb(ref); ms=pdbqt_models(dock)
    except Exception: continue
    if len(cx)==0 or len(ms)<2: continue
    cases.append((cid,cx,ms))
print(f"usable cases: {len(cases)}")

# crystal vectors
crystal={cid:pv.contact_vector(cx,rec,RESI) for cid,cx,_ in cases}

hits=0; tot=0; ranks=[]
for cid,cx,ms in cases:
    # LEAVE-ONE-OUT reference profile: this case's own crystal pose excluded
    ref_vs=[v for k,v in crystal.items() if k!=cid]
    med,mad=pv.reference_profile(ref_vs, threshold=3.0)
    cand=[crystal[cid]]+[pv.contact_vector(m,rec,RESI) for m in ms]
    scores=[pv.fit_score(v,med,mad) for v in cand]
    order=np.argsort(scores)
    rank=int(np.where(order==0)[0][0])+1     # rank of the CRYSTAL pose
    ranks.append(rank); tot+=1
    if rank==1: hits+=1
print()
print(f"crystal pose ranked #1 of {len(cand)} by fit score: {hits}/{tot} = {100*hits/tot:.1f}%")
print(f"median rank of the crystal pose: {np.median(ranks):.1f}  (chance = {(len(cand)+1)/2:.1f})")
print(f"crystal in top 3: {sum(1 for r in ranks if r<=3)}/{tot} = {100*sum(1 for r in ranks if r<=3)/tot:.1f}%")
