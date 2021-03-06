import numpy as np
import sys
sys.path.insert(0, '/home/david/Work/DFT-Surface-Science/sublattice-method/funcs')
from funcs import rep, rot, lsave, genPOSCAR

lvecs = np.genfromtxt('../in/_lvecs')
bases=np.genfromtxt('../in/_bases')
zn1,zn2,o1,o2=bases
na=1
nb=1
nc=1

lstzn1=rep(zn1,lvecs,na,nb,nc)
lstzn2=rep(zn2,lvecs,na,nb,nc)
lsto1=rep(o1+lvecs[0],lvecs,na,nb,nc)
lsto2=rep(o2,lvecs,na,nb,nc)

lsave(lstzn1,'Zn1')
lsave(lstzn2,'Zn2')
lsave(lsto1,'O1')
lsave(lsto2,'O2')

sublats=['Zn1','Zn2','O1','O2']

a=np.array([lvecs[0],np.zeros(3),np.zeros(3)])
b=np.array([np.zeros(3),lvecs[1],np.zeros(3)])
c=np.array([np.zeros(3),np.zeros(3),lvecs[2]])
rng=np.arange(0.85,1.15,0.0125)
c=0
for i in rng:
    for j in rng:
            c+=1
            scvecs=[lvecs[0]*i,lvecs[1],lvecs[2]*j]
            genPOSCAR(sublats,scvecs,sys_name='nonpolar zinc oxide',POS_name=str(c)+'acPOSCAR.vasp')
