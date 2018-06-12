import numpy as np
import sys
sys.path.insert(0, '/home/david/Work/DFT-Surface-Science/sublattice-method/funcs')
from funcs import orient, rotation_matrix, lsave, genPOSCAR
import csv
H_lst=[]
O_lst=[]
N_lst=[]
C_lst=[]

with open('../in/_bases','r') as csvfile:
    csvreader= csv.reader(csvfile,delimiter=' ',skipinitialspace=True)
    for line in csvreader:
        x,y,z,name,_=line
        x=float(x)
        y=float(y)
        z=float(z)
        dist=np.array([x,y,z])
        if name=='H':
            H_lst.append(dist)
        elif name=='O':
            O_lst.append(dist)
        elif name=='N':
            N_lst.append(dist)
        else:
            C_lst.append(dist)
H_rot=[]
O_rot=[]
N_rot=[]
C_rot=[]
#rotation is not about the center of the molecule but about the origin, need to translate it after the rotation
#the molecule will be incorrectly rotated if it is not reported in terms of coordinates relative to its. This may be found
alt=90*np.pi/180
azi=90*np.pi/180
disp=np.array([5,5,5])

H_lst=orient(H_lst,alt,azi)+disp
O_lst=orient(O_lst,alt,azi)+disp
N_lst=orient(N_lst,alt,azi)+disp
C_lst=orient(C_lst,alt,azi)+disp

lsave(H_lst,'H1')
lsave(O_lst,'O1')
lsave(N_lst,'N1')
lsave(C_lst,'C1')

sublats=['H1','O1','N1','C1']

scvecs=[[10,0,0],[0,10,0],[0,0,10]]

genPOSCAR(sublats,scvecs,sys_name='molecule',POS_name='mculPOSCAR.vasp')
