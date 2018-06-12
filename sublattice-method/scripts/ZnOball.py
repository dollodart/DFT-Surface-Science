import numpy as np
import sys
sys.path.insert(0, '/home/david/Work/DFT-Surface-Science/sublattice-method/funcs')
from funcs import rep, rot, lsave, genPOSCAR

lvecs = np.genfromtxt('../in/_lvecs')
bases=np.genfromtxt('../in/_bases')
zn1,zn2,o1,o2=bases
a,b,c=lvecs
R=12. #Angstroms
arange=int(np.ceil(R/np.linalg.norm(a)))
brange=int(np.ceil(R/np.linalg.norm(b)))
crange=int(np.ceil(R/np.linalg.norm(c)))
comb=[]
lstzn1=[]
lstzn2=[]
lsto1=[]
lsto2=[]

i1=0
i2=0
for na in range(-arange,arange):
    for nb in range(-brange,brange):
        for nc in range(-crange,crange):
            dist = np.linalg.norm(na*a+nb*b+nc*c)
            i1+=1
            if dist < R:
                i2+=1
                comb.append([na,nb,nc])

for i in comb:
    dist=np.dot(np.array(i),lvecs)+np.array([2*R,2*R,2*R])
    lstzn1.append(zn1+dist)
    lstzn2.append(zn2+dist)
    lsto1.append(o1+a+dist)
    lsto2.append(o2+dist) #since don't account for non-centro symmetry, get non-physical terminations. Would need logical similar to that of the slab replication

lsave(lstzn1,'Zn1') 
lsave(lstzn2,'Zn2')
lsave(lsto1,'O1')
lsave(lsto2,'O2')

sublats=['Zn1','Zn2','O1','O2']

scvecs=[[4*R,0,0],[0,4*R,0],[0,0,4*R]]
genPOSCAR(sublats,scvecs,sys_name='zinc oxide ball',POS_name='ZnOballPOSCAR.vasp')
