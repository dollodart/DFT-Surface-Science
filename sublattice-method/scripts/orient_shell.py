import numpy as np
import sys
sys.path.insert(0, '/home/david/Work/DFT-Surface-Science/sublattice-method/funcs')
from funcs import *

R=10. #Angstroms

mculo2=np.array([[0,0,0],[-1.208,0,0]]) #diatomic oxygen. Easier than a molecule which is chiral or otherwise has not only an orientation but a head and a tail

unit_x=np.array([1.,0,0])
unit_y=np.array([0,1.,0])
unit_z=np.array([0,0,1.])

a=5.*unit_x
b=5.*unit_y
c=5.*unit_z
lvecs=np.array([a,b,c])

arange=int(np.ceil(R/np.linalg.norm(a)))
brange=int(np.ceil(R/np.linalg.norm(b)))
crange=int(np.ceil(R/np.linalg.norm(c)))

shell=[]

for na in range(-arange,arange+1):#take into account up-to-and-excluding loop convention
    for nb in range(-brange,brange+1):
        for nc in range(-crange,crange+1):
            vct=na*a+nb*b+nc*c
            r=np.linalg.norm(vct)
            if r < R+2 and r > R-2: 
                shell.append([na,nb,nc]) 

lstmcul=[]
epsilon=1.e-3
for i in shell:
    vct=np.dot(np.array(i),lvecs)
    r=np.linalg.norm(vct)
    # only sign ambiguities at this point
    if abs(np.dot(vct,unit_z)) < epsilon:
        alt=np.pi/2
    else:
        alt = np.arctan(np.linalg.norm(vct[0:2])/np.dot(vct,unit_z))
    if abs(np.dot(vct,unit_x)) < epsilon:
        azi=np.pi/2
    else:
        azi = np.arctan(np.dot(vct,unit_y)/np.dot(vct,unit_x))
    new=orient(mculo2,azi,alt)
    for atom in new:
        tot=atom+vct+np.array([1,1,1])*2*R
        lstmcul.append(tot) 
lsave(lstmcul,'O1')

sublats=['O1']

scvecs=[[4*R,0,0],[0,4*R,0],[0,0,4*R]]

genPOSCAR(sublats,scvecs,sys_name='adsorbate oriented shell',POS_name='orientshellPOSCAR.vasp')
