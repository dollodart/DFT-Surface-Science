import numpy as np
def rep(bas,vecs,na,nb,nc):
    """replicates the bases vector by lattice vectors to the specified multiples"""
    acc=[]
    for i in range(na):
	    for j in range(nb):
		    for k in range(nc):
			    mult=np.array([i,j,k])
			    pos=np.dot(mult,vecs)+bas			
			    acc.append(pos)
    return acc

#def rot(vecs,z,y,x): # about z, about y, about x
#    z=z*np.pi/180
#    y=y*np.pi/180
#    x=x*np.pi/180
#
#    cx=np.cos(x)
#    sx=np.sin(x)
#    cy=np.cos(y)
#    sy=np.sin(y)
#    cz=np.cos(z)
#    sz=np.sin(z)
#
#    Rx=np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]])
#    Ry=np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]])
#    Rz=np.array([[cz,sz,0],[-sz,cz,0],[0,0,1]])
#    R=np.matmul(Rx,np.matmul(Ry,Rz))
#    rvecs=[]
#
#    for vec in vecs:
#            rvec=np.dot(vec,R) #left multiplying row vector		
#            rvecs.append(rvec)
#    return rvecs
# above only works for orthogonal, independent angular transforms. 

def rot(vecs,z,y,x): # about z, about y, about x
    z=z*np.pi/180
    y=y*np.pi/180
    x=x*np.pi/180
    c1=np.cos(z)
    c2=np.cos(y)
    c3=np.cos(x)
    s1=np.sin(z)
    s2=np.sin(y)
    s3=np.sin(x)
    r11=c1*c2
    r12=c1*s2*s3-c3*s1
    r13=s1*s3+c1*c3*s2
    r21=c2*s1
    r22=c1*c3+s1*s2*s3
    r23=c3*s1*s2-c1*s3
    r31=-s2
    r32=c2*s3
    r33=c2*c3
    r=[[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]]
    rvecs=[]
    for vec in vecs:
        rvec=np.dot(vec,r)
        rvecs.append(rvec)
    return rvecs

ddir='../.dump/'

def lsave(vecs,fname):
    """ saves the lattice to storage """
    np.savetxt(ddir+fname,vecs)
def genPOSCAR(sublats,scvecs,sys_name=None,POS_name=None):
    """ Generates the poscar file """
    sys_name = sys_name or 'no name given'
    POS_name = POS_name or 'POSCAR'
    lat=[]
    nums=[]
    specs=[]
    hist=[]

    tt=dict.fromkeys(map(ord, '1234567890'),None)

    for sublat in sublats:
        sl=np.genfromtxt(ddir+sublat)
        sln=int(sl.size/3) # 3 coordinates per atom
        lat.append(sl)
        spec=sublat.translate(tt)
        if spec in hist:
            ari=float(nums[-1])+sln
            nums[-1]=str(int(ari))+' '
            # do not append a new species
        else:
            specs.append(spec)
            nums.append(str(sln)+' ')
        hist.append(spec)

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

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def orient(pos0s,azi,alt):
    """ Return the orientation outward normal with angles specified in radians"""
    axis1 = np.array([0,0,1]) # z-axis
    axis2 = np.array([-np.sin(azi),np.cos(azi),0]) # vector perpendicular to projection in x-y plane. Note orientation only one orientation of the two possible for the axis is selected, it being a vectory quantity, and so it may fail symmetry

    if len(pos0) > 1:
        pos=[]
        for pos0 in pos0s:
            pos1 = np.dot(rotation_matrix(axis1, azi),pos0)
            pos2 = np.dot(rotation_matrix(axis2, alt),pos1)
            pos.append(pos2)
        return pos

    else:
        pos1 = np.dot(rotation_matrix(axis1, azi),pos0)
        pos2 = np.dot(rotation_matrix(axis2, alt),pos1)
        pos.append(pos2)
