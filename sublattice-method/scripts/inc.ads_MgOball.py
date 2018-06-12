import numpy as np
import sys
sys.path.insert(0, '/home/david/Work/DFT-Surface-Science/sublattice-method/funcs')
from funcs import *

R=15. #Angstroms
lvecs = np.genfromtxt('../in/_lvecs')
a,b,c=lvecs

bases=np.genfromtxt('../in/_bases')
mg1,mg2,mg3,mg4,o1,o2,o3,o4=bases
mculo2=[[0,0,0],[-1.208,0,0]] #diatomic oxygen

arange=int(np.ceil(R/np.linalg.norm(a)))
brange=int(np.ceil(R/np.linalg.norm(b)))
crange=int(np.ceil(R/np.linalg.norm(c)))
comb=[]
shell=[]

lstmg1=[]
lstmg2=[]
lstmg3=[]
lstmg4=[]
lsto1=[]
lsto2=[]
lsto3=[]
lsto4=[]

ncube=8*arange*brange*crange

nsph=0
i=np.array([1,0,0])
j=np.array([0,1,0])
k=np.array([0,0,1])


for na in range(-arange,arange):
    for nb in range(-brange,brange):
        for nc in range(-crange,crange):
            vct=na*a+nb*b+nc*c
            dist = np.linalg.norm(vct)
            if dist < R:
                comb.append([na,nb,nc])
                nsph+=1
            if dist > R and dist < R + 2: # need to know increment of length
                shell.append([na,nb,nc]) 

for i in comb: #figure out how to initialize a 2-D structure
    dist=np.dot(np.array(i),lvecs) + np.array([2*R,2*R,2*R])
    lstmg1.append(mg1+dist)
    lstmg2.append(mg2+dist)
    lstmg3.append(mg3+dist)
    lstmg4.append(mg4+dist)
    lsto1.append(o1+dist)
    lsto2.append(o2+dist)
    lsto3.append(o3+dist)
    lsto4.append(o4+dist)
lstmcul=[]
for molecule in shell:
    vct=np.dot(np.array(i),lvecs)
    dist=np.linalg.norm(vct)
    alt=np.arccos(np.dot(vct,k)/dist)
    azi=np.arctan(np.dot(vct,j)/np.dot(vct,i))*180/np.pi
    for atom in molecule:
        pos=orient(atom,alt,azi)
        lstmcul.append(pos) #mcul is made of many atoms and so will require better list management

lsave(lstmg1,'Mg1')  #array method possible?
lsave(lstmg2,'Mg2')
lsave(lstmg3,'Mg3') 
lsave(lstmg4,'Mg4')
lsave(lsto1,'O1')
lsave(lsto2,'O2')
lsave(lsto3,'O3')
lsave(lsto4,'O4')
lsave(lstmcul,'P1') # for visual distinction

sublats=['Mg1','Mg2','Mg3','Mg4','O1','O2','O3','O4','P1']

scvecs=[[4*R,0,0],[0,4*R,0],[0,0,4*R]]

genPOSCAR(sublats,scvecs,sys_name='Adsorbate Magnesium Oxide Ball System',POS_name='ads_MgOBallPOSCAR.vasp')
metr=100*(nsph/ncube-np.pi/6)/(np.pi/6)
print(metr)
