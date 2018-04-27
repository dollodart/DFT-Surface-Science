import numpy as np
import sys
sys.path.insert(0, '/root/temp/funcs') # change to pwd of machine
from funcs import rep, rot, lsave, genPOSCAR
# eventually fix array dimensioning

lvecs = np.genfromtxt('in/_lvecs')

bases=np.genfromtxt('in/_bases')
lst=[]
for i in range(8):
	lst.append(bases[i])

lst2=[]
for i in lst:
	lst2.append(rep(i,lvecs,1,1,5))
# rot(lst,yaw,pitch,roll) will rotate vectors about the coordinate system

sublats=['Mg1', 'Mg2', 'Mg3', 'Mg4', 'O1', 'O2', 'O3', 'O4']
o=0
for i in sublats:
	lsave(lst2[o],i)
	o+=1

scvecs=lvecs + np.array([[0,0,0],[0,0,0],[0,0,15+lvecs[2,2]*4]]) # vacuum spacing

genPOSCAR(sublats,scvecs,'Magnesium Oxide 4 Slab System')
