import numpy as np
import sys
sys.path.insert(0, '/home/david/Work/DFT-Surface-Science/sublattice-method/funcs')
from funcs import rep, rot, lsave, genPOSCAR

lvecs = np.genfromtxt('../in/_lvecs')
bases=np.genfromtxt('../in/_bases')
zn1,zn2,o1,o2=bases
na=1
nb=1

for nc in range(3,11):
    ncl=int(np.floor(nc/2))
    ncu=int(np.ceil(nc/2))
    shift=nc/4.*lvecs[2]

    if nc % 2 == 0: #even
        lstzn1=rep(zn1+shift,lvecs,na,nb,ncu) 
        lstzn2=rep(zn2+shift,lvecs,na,nb,ncl-1)
        lsto1=rep(o1+lvecs[2]+shift,lvecs,na,nb,ncl-1)
        lsto2=rep(o2+shift,lvecs,na,nb,ncl)
    else: # odd
        lstzn1=rep(zn1+shift,lvecs,na,nb,ncl)
        lstzn2=rep(zn2+shift,lvecs,na,nb,ncl)
        lsto1=rep(o1+lvecs[2]+shift,lvecs,na,nb,ncl)
        lsto2=rep(o2+shift,lvecs,na,nb,ncl)

    lsave(lstzn1,'Zn1')
    lsave(lstzn2,'Zn2')
    lsave(lsto1,'O1')
    lsave(lsto2,'O2')

    Zn_Ot=o1+(nc/2-0.3)*lvecs[2] + shift
    O_Ot=zn2+(nc/2-1.1)*lvecs[2] + shift
    hole_Ot=lvecs[0]+lvecs[1]+(nc/2-0.3)*lvecs[2] + shift
    bridge_Ot=lvecs[0]/2+lvecs[1]/2+(nc/2-0.3)*lvecs[2] + shift 

    Zn_Znt=zn2-0.3*lvecs[2] + shift 
    O_Znt=zn1-0.3*lvecs[2] + shift
    hole_Znt=lvecs[0]+lvecs[1]+0.3*lvecs[2] + shift
    bridge_Znt=lvecs[0]/2+lvecs[1]/2+0.3*lvecs[2] + shift

    lstO=rep(O_Znt-0.1*lvecs[2],lvecs,na,nb,1)
    lstH_Znt=rep(O_Znt-0.325*lvecs[2],lvecs,na,nb,1) 
    if nc % 2 == 0:
        lstH_Ot=rep(O_Ot-0.0625*lvecs[2],lvecs,na,nb,1)
    else:
        lstH_Ot=rep(Zn_Ot,lvecs,na,nb,1)
    lsave(lstH_Znt,'H1')
    lsave(lstO,'O3')
    lsave(lstH_Ot,'H2')

    sublats=['Zn1','Zn2','O1','O2','O3','H1','H2']

    scvecs=[na*lvecs[0],nb*lvecs[1],nc*lvecs[2]]

    genPOSCAR(sublats,scvecs,sys_name='Zinc Oxide H-terminated',POS_name=str(nc)+'HpPOSCAR.vasp')
