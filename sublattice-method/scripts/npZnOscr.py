import numpy as np
import sys
sys.path.insert(0, '/root/temp/funcs') # change to pwd of machine
from funcs import rep, rot, lsave, genPOSCAR

lvecs = np.genfromtxt('../in/_lvecs')

bases=np.genfromtxt('../in/_bases')

zn1,zn2,o1,o2=bases

na=10
nb=2
nc=2 

lstzn1=rep(zn1,lvecs,na,nb,nc)
lstzn2=rep(zn2,lvecs,na,nb,nc)
lsto1=rep(o1+lvecs[0],lvecs,na,nb,nc)
lsto2=rep(o2,lvecs,na,nb,nc)

#note that rot(lst,yaw,pitch,roll) will rotate vectors about the coordinate system

lsave(lstzn1,'Zn1') #convention is to use numerical end characters with the proper species name prefix
lsave(lstzn2,'Zn2')
lsave(lsto1,'O1')
lsave(lsto2,'O2')

sublats=['Zn1','Zn2','O1','O2']

a=np.array([lvecs[0],np.zeros(3),np.zeros(3)])
b=np.array([np.zeros(3),lvecs[1],np.zeros(3)])
c=np.array([np.zeros(3),np.zeros(3),lvecs[2]])
scvecs=na*a+nb*b+nc*c+4*a
genPOSCAR(sublats,scvecs,sys_name='nonpolar zinc oxide',POS_name='npPOSCAR.vasp')
