import numpy as np
import sys
sys.path.insert(0, '/home/david/Work/DFT-Surface-Science/sublattice-method/funcs')
from funcs import rep, rot, lsave, genPOSCAR

lvecs = np.genfromtxt('../in/_lvecs')
bases=np.genfromtxt('../in/_bases')
zn1,zn2,o1,o2=bases
nar = range(3,20)

for na in nar:
    nb=2
    nc=2 

    lstzn1=rep(zn1,lvecs,na,nb,nc)
    lstzn2=rep(zn2,lvecs,na,nb,nc)
    lsto1=rep(o1+lvecs[0],lvecs,na,nb,nc)
    lsto2=rep(o2,lvecs,na,nb,nc)

    lsave(lstzn1,'Zn1') #convention is to use numerical end characters with the proper species name prefix
    lsave(lstzn2,'Zn2')
    lsave(lsto1,'O1')
    lsave(lsto2,'O2')

    O2=o2+(na-1/3)*lvecs[0]
    Zn2=zn2+(na-1/3)*lvecs[0]
    O1=o1+(na+1/3)*lvecs[0]
    Zn1=zn1+(na-1/3)*lvecs[0]
    hole=(o1+zn1)/2+1/2*lvecs[1] + na*lvecs[0]
    bridge=o2+(na-1/2)*lvecs[0]+lvecs[2]/5
    hole2=zn2+(lvecs[2]+lvecs[1])/4 + na*lvecs[0] #off center, subtract c to get at same position as hole

    lsave(O2,'S1')
    lsave(Zn2,'S2')
    lsave(O1,'S3')
    lsave(Zn1,'S4')
    lsave(hole,'S5')
    lsave(bridge,'S6')
    lsave(hole2, 'S7')

    sublats=['Zn1','Zn2','O1','O2','S1','S2','S3','S4','S5','S6','S7']

    scvecs=[(na+4)*lvecs[0],nb*lvecs[1],nc*lvecs[2]]
    genPOSCAR(sublats,scvecs,sys_name='nonpolar zinc oxide',POS_name=str(na)+'npPOSCAR.vasp')
