import numpy as np
def rep(bas,vecs,na,nb,nc):
	acc=[]
	for i in range(na):
		for j in range(nb):
			for k in range(nc):
				mult=np.array([i,j,k])
				pos=np.dot(mult,vecs)+bas			
				acc.append(pos)
	return acc

def rot(vecs,yaw,pitch,roll): # about z, about y, about x
	yaw=yaw*np.pi/180
	pitch=pitch*np.pi/180
	roll=roll*np.pi/180

	cx=np.cos(roll)
	sx=np.sin(roll)
	cy=np.cos(pitch)
	sy=np.sin(pitch)
	cz=np.cos(yaw)
	sz=np.sin(yaw)

	Rx=np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]])
	Ry=np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]])
	Rz=np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]])
	R=np.matmul(Rx,np.matmul(Ry,Rz))
	rvecs=[]

	for vec in vecs:
		rvec=np.dot(vec,R) #left multiplying row vector		
		rvecs.append(rvec)
	return rvecs
ddir='../.dump/'
def lsave(vecs,fname):
	np.savetxt(ddir+fname,vecs)
def genPOSCAR(sublats,scvecs,sys_name=None,POS_name=None):
	sys_name = sys_name or 'no name given'
	POS_name = POS_name or 'POSCAR'
	lat=[]
	nums=[]
	specs=[]

	tt=dict.fromkeys(map(ord, '1234567890'),None)

	for sublat in sublats:
		sl=np.genfromtxt(ddir+sublat)
		sln=int(sl.size/3) # 3 coordinates per atom
		lat.append(sl)
		nums.append(str(sln)+' ')
		specs.append(sublat.translate(tt)+' ')
	stack=np.vstack(lat)
	np.savetxt(ddir + 'lat',stack)
	np.savetxt(ddir + 'scvecs',scvecs)

	g=open(ddir + 'lat','r')
	h=open(ddir + 'scvecs','r')
	f=open('../out/' + POS_name,'w')
	f.write(sys_name + '\n')
	f.write('1.0\n')

	for line in h: 
		f.write(line)
	for spec in specs:
		f.write(spec + ' ')
	f.write('\n')
	for num in nums:
		f.write(num + ' ')
	f.write('\nCartesian\n')
	for line in g:
		f.write(line)
	h.close()
	g.close()
	f.close()
