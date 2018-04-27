import numpy as np
import math

vecs = np.array([[3.2891030312,0,0],[-1.6445515156,2.8484467806,0],[0,0,5.3068203926]])

zn1=np.array([1.644553185,0.949481308,2.650501966])
zn2=np.array([-0.0,1.898965478,5.303912163])
o1=np.array([-1.644553185,0.949481308,0.638081491])
o2=np.array([-0.0,1.898965478,3.291491747])

def sl(bas,fname,mm,nn,p):
	ap=[]
	for n in range(nn):
		for m in range(mm):
#			rndpert=np.random.rand(3)/100 # visualization sensitive to perturbations 
			rndpert=np.zeros(3)
			ap.append(bas+np.dot(np.array([m,n,p/2]),vecs) + rndpert)
	np.savetxt(fname,np.squeeze(ap))
def nsl(i,f):
	g = open(f,'w')
	g.write(str(i))
	g.close()	

mm=2
nn=2
pp=20

dr='pslabs/d'

tbas=o1+np.array([0,0,1.5])
bbas=zn1+np.array([0,0,-1.5])
for p in range(pp):
	sl(tbas,dr + str(p) + 'tads',mm,nn,p)
	sl(bbas,dr + str(p) + 'bads',mm,nn,0)
	nsl(mm*nn,dr + str(p) + 'tadsnum') 
	nsl(mm*nn,dr + str(p) + 'badsnum') 


nsl('S S','pslabs/ads_id')
