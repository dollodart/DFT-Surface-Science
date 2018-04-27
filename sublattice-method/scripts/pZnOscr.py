import numpy as np
import sys
sys.path.insert(0, '/root/temp/funcs') # change to pwd of machine
from funcs import rep, rot, lsave, genPOSCAR

lvecs = np.genfromtxt('../in/_lvecs')

bases=np.genfromtxt('../in/_bases')

zn1,zn2,o1,o2=bases

na=2
nb=2
nc=7 

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

#note that rot(lst,yaw,pitch,roll) will rotate vectors about the coordinate system

lsave(lstzn1,'Zn1') #convention is to use numerical end characters with the proper species name prefix
lsave(lstzn2,'Zn2')
lsave(lsto1,'O1')
lsave(lsto2,'O2')

sublats=['Zn1','Zn2','O1','O2']

a=np.array([lvecs[0],np.zeros(3),np.zeros(3)])
b=np.array([np.zeros(3),lvecs[1],np.zeros(3)])
c=np.array([np.zeros(3),np.zeros(3),lvecs[2]])
scvecs=na*a+nb*b+nc*c
genPOSCAR(sublats,scvecs,sys_name='Zinc Oxide',POS_name='pPOSCAR.vasp')
