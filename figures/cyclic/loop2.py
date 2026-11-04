import sys,re,numpy as np
from glob import glob
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent.parent/"port_2611"))
from ensight_reader_2611 import read_geo, read_var_node_all_steps
run=Path(sys.argv[1]); h=sys.argv[2]; angle=10
d=run/(("Cohesive_zone_model" if (run/f"Cohesive_zone_model_top_{h}").exists() else "Cauchy_Full_model")+f"_top_{h}")
parts=read_geo(str(d/"geometry.geo_0000"))
def arr(var):
    fs=sorted(glob(str(d/f"{var}_top.var_*")));g=[read_var_node_all_steps(f,parts,3) for f in fs];pn=list(g[0])[0];return np.array([x[pn][-1] for x in g])
R=arr("reaction");U=arr("displacement")
case=Path(str(d)+".case").read_text()
t=np.array(re.findall(r"[-+0-9.e]+",case.split("time values:")[-1].split("GEOMETRY")[0])[:len(R)],float)
a=np.radians(angle);dl=np.array([np.cos(a),0,-np.sin(a)])
F=-np.einsum('tni,i->t',R,dl); u=np.einsum('tni,i->tn',U,dl).mean(axis=1)
np.savez(run/"loop.npz",t=t,u=u,F=F)
for ti,uu,f in zip(t,u,F): print(f"t={ti:6.3f} u={uu: .3e} F={f: .4e}")
