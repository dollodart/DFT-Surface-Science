import numpy as np
import sys
sys.path.insert(0, '/home/david/Work/DFT-Surface-Science/sublattice-method/funcs')
from funcs import rep, rot, lsave, genPOSCAR
ncr=range(3,20)

for nc in ncr:
    lvecs = np.genfromtxt('../in/_lvecs')
    bases=np.genfromtxt('../in/_bases')
    zn1,zn2,o1,o2=bases
    na=4
    nb=4
    ncl=int(np.floor(nc/2))
    ncu=int(np.ceil(nc/2))

    if nc % 2 == 0: #even
            lstzn1=rep(zn1,lvecs,na,nb,ncu)
            lstzn2=rep(zn2,lvecs,na,nb,ncl-1)
            lsto1=rep(o1+lvecs[2],lvecs,na,nb,ncl-1)
            lsto2=rep(o2,lvecs,na,nb,ncl)
    else: # odd
            lstzn1=rep(zn1,lvecs,na,nb,ncl)
            lstzn2=rep(zn2,lvecs,na,nb,ncl)
            lsto1=rep(o1+lvecs[2],lvecs,na,nb,ncl)
            lsto2=rep(o2,lvecs,na,nb,ncl)

    lsave(lstzn1,'Zn1')
    lsave(lstzn2,'Zn2')
    lsave(lsto1,'O1')
    lsave(lsto2,'O2')
    #bases of surface sites
    #have arbitrary positions above surfaces, that is, multiples of lvecs[2]
    #due to non-centro symmetry, capping will change on even and odd numbered layers. These names are given for even number of double layers, opposite for odd number of double layers
    Zn_Ot=o1+(nc/2-0.3)*lvecs[2]
    O_Ot=zn2+(nc/2-1.1)*lvecs[2]
    hole_Ot=lvecs[0]+lvecs[1]+(nc/2-0.3)*lvecs[2] 
    bridge_Ot=lvecs[0]/2+lvecs[1]/2+(nc/2-0.3)*lvecs[2] 

    Zn_Znt=zn2-0.3*lvecs[2] 
    O_Znt=zn1-0.3*lvecs[2] 
    hole_Znt=lvecs[0]+lvecs[1]+0.3*lvecs[2]
    bridge_Znt=lvecs[0]/2+lvecs[1]/2+0.3*lvecs[2]

    lsave(hole_Ot,'S1')
    lsave(bridge_Ot,'S2')
    lsave(Zn_Ot,'S4')
    lsave(O_Ot,'S3')

    lsave(hole_Znt,'S5')
    lsave(bridge_Znt,'S6')
    lsave(O_Znt,'S7')
    lsave(Zn_Znt,'S8')

    sublats=['Zn1','Zn2','O1','O2','S1','S2','S3','S4','S5','S6','S7','S8']

    scvecs=[na*lvecs[0],nb*lvecs[1],nc*lvecs[2]]
    genPOSCAR(sublats,scvecs,sys_name='Zinc Oxide',POS_name=str(nc)+'pPOSCAR.vasp')
