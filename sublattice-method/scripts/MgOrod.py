import numpy as np
import sys
sys.path.insert(0, '/home/david/Work/DFT-Surface-Science/sublattice-method/funcs')
from funcs import rep, rot, lsave, genPOSCAR

R=20. #Angstroms
x=2
y=2
z=8
lvecs = np.genfromtxt('../in/_lvecs')
a,b,c=lvecs
bases=np.genfromtxt('../in/_bases')
arange=int(np.ceil(R/np.linalg.norm(a)))
brange=int(np.ceil(R/np.linalg.norm(b)))
crange=8
mg1,mg2,mg3,mg4,o1,o2,o3,o4=bases
comb=[]
lstmg1=[]
lstmg2=[]
lstmg3=[]
lstmg4=[]
lsto1=[]
lsto2=[]
lsto3=[]
lsto4=[]
nrect=8*arange*brange*crange
ncyl=0
for na in range(-arange,arange):
    for nb in range(-brange,brange):
        for nc in range(-crange,crange):
            dist = np.linalg.norm(na*a+nb*b)
            if dist < R:
                comb.append([na,nb,nc])
                ncyl+=1

for i in comb:
    dist=np.dot(np.array(i),lvecs) + np.array([2*R,2*R,0])
    lstmg1.append(mg1+dist)
    lstmg2.append(mg2+dist)
    lstmg3.append(mg3+dist)
    lstmg4.append(mg4+dist)
    lsto1.append(o1+dist)
    lsto2.append(o2+dist)
    lsto3.append(o3+dist)
    lsto4.append(o4+dist)

lsave(lstmg1,'Mg1')  #array method possible?
lsave(lstmg2,'Mg2')
lsave(lstmg3,'Mg3') 
lsave(lstmg4,'Mg4')
lsave(lsto1,'O1')
lsave(lsto2,'O2')
lsave(lsto3,'O3')
lsave(lsto4,'O4')

sublats=['Mg1','Mg2','Mg3','Mg4','O1','O2','O3','O4']

scvecs=[[4*R,0,0],[0,4*R,0],[0,0,(2*crange+4)*np.linalg.norm(c)]]

genPOSCAR(sublats,scvecs,sys_name='Magnesium Oxide Rod System',POS_name='MgORodPOSCAR.vasp')
rel=np.pi*crange*np.linalg.norm(c)/R
metr=100*(ncyl/nrect-rel)/rel
print(metr)
