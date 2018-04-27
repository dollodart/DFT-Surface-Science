import numpy as np
import math

vecs = np.array([[0,0,5.3068203926],[-1.6445515156,2.8484467806,0],[3.2891030312,0,0]])

def sl(bas,fname,mm,nn,pp,pp0):
	ap=[]
	for p in range(pp0,pp):
		for n in range(nn):
			for m in range(mm):
#				rndpert=np.random.rand(3)/100 # visualization in VESTA can be sensitive to small perturbations
				rndpert=np.zeros(3)
				ap.append(bas+np.dot(np.array([m,n,p]),vecs) + rndpert)
	np.savetxt(fname,np.squeeze(ap))
def nsl(i,f):
	g = open(f,'w')
	g.write(str(i))
	g.close()	


zn1=np.array([1.644553185,0.949481308,2.650501966])
zn2=np.array([-0.0000,1.898965478,5.303912163])
o1=np.array([-1.644553185,0.949481308,0.638081491])
o2=np.array([-0.0000,1.898965478,3.291491747])

#dimensions of supercell, will generate successive slabs to pp-1. The shell script must be updated for larger pp
mm=2
nn=2
pp=20

avecs=vecs.tolist()
bvecs=np.array(avecs)

dr='npslabs/d'

bvecs[0,:]=bvecs[0,:]*mm
bvecs[1,:]=bvecs[1,:]*nn
bvecs[2,0]=bvecs[2,0]+15 #vacuum spacing of 15 Angstroms 

for p in range(pp):
	ps=p+1 
	pe=ps/2
	phigh=math.ceil(pe)
	plow=math.floor(pe)
	if ps % 2 == 0:
		sl(zn1,dr + str(p) + 'zn1',mm,nn,phigh,0)
		sl(zn2,dr + str(p) + 'zn2',mm,nn,plow-1,0)
		sl(o1,dr + str(p) + 'o1',mm,nn,phigh,1)
		sl(o2,dr + str(p) + 'o2',mm,nn,phigh,0)
	if ps % 2 == 1:
		sl(zn1,dr + str(p) + 'zn1',mm,nn,plow,0)
		sl(zn2,dr + str(p) + 'zn2',mm,nn,phigh-1,0)
		sl(o1,dr + str(p) + 'o1',mm,nn,phigh,1)
		sl(o2,dr + str(p) + 'o2',mm,nn,plow,0)
		
	nsl((ps-1)*mm*nn,dr + str(p) + 'num') #subtraction to account for termination eliminations 
	bvecs[2,0]=bvecs[2,0] + vecs[2,0]/2
	np.savetxt(dr + str(p) + 'vecs',bvecs)
