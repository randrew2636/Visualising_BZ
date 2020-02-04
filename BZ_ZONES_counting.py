# Richard Andrew 2015
# 
import numpy as np
import pandas as pd

import time
import os
#--------------------------------------------------------

if (not os.path.exists('ws.cpython-37m-darwin.so')):
    os.system("f2py3.7 -c -m ws bzones_lib_counting.f90 --f90flags=-fopenmp -lgomp")

import ws

# Rvect generator
def gvect(gs,gcut,n):
    # define spanning vectors
    g = np.zeros(((n+2)*(n+2)*(n+2)*8*8*8,3))
    l=0
    # Generate all G vectors in more or less order of increasing shell size within cutoff radius
    print ('Generating Gvect')
    for i in range(-n,n+1):
        for j in range(-n,n+1):
            for k in range(-n,n+1):
                if(i==0 and j==0 and k==0):
                    continue
                g[l,:] = i*gs[0,:] + j*gs[1,:] + k*gs[2,:]
                gg = np.linalg.norm(g[l])**2
                if(gg-gcut > 0.0001):
                    continue
                l=l + 1
    print ()
    print ('Number of vectors= ',l)
    # order G vectors into shells of increasing radius
    gg = np.zeros(l)
    for i in range(l):
        gg[i]=np.linalg.norm(g[i])**2
    idx = np.argsort(gg)
    gg=gg[idx]
    g=g[idx]
    g=g[0:l]
    return g
#----------------------------------------------------------------------

PI = np.pi
#### Enter in parameters
# Enter crystal system type
print ('Enter cystal system type:')
print ()
print ('1-CUBIC')
print ('2-TETRAGONAL')
print ('3-ORTHORHOMBIC')
print ('4-MONOCLINIC')
print ('5-TRICLINIC')
print ('6-HEXAGONAL')
print ('7-TRIGONAL')
print ()
crys = int(input('Enter: '))

# Enter crystal type and associated parameters
ca=1.0
ba=1.0
if(crys==1):
    print ('P-1 I-2 F-3')
    s = int(input('Enter: '))
elif(crys==2):
    print ('P-1 I-2')
    s = int(input('Enter: '))
    ca = float(input('enter c/a '))
elif(crys==3):
    print ('P-1 I-2 F-3 C-4')
    s = int(input('Enter: '))
    ba = float(input('enter b/a '))
    ca = float(input('enter c/a '))
elif(crys==4):
    print ('P-1 I-2-NOT CHECKED')
    s = int(input('Enter: '))
    ba = float(input('enter b/a '))
    ca = float(input('enter c/a '))
    gamma = float(input('enter gamma angle between a anc b vectors '))
    gamma=gamma*PI/180.0
elif(crys==5):
    ba = float(input('enter b/a '))
    ca = float(input('enter c/a '))

    alpha = float(input('enter alpha angle between b anc c vectors '))
    beta  = float(input('enter beta  angle between a anc c vectors '))
    gamma = float(input('enter gamma angle between a anc b vectors '))

    alpha=alpha*PI/180.0
    beta=beta*PI/180.0
    gamma=gamma*PI/180.0
else:
    ca = float(input('enter c/a '))

print ()

# RS are the spanning vectors
# R are the generated lattice vectors
RS = np.zeros(shape=(3,3))


# Set the spanning vectors in units of 'a'
if (crys==1 and s==1):
    RS[0,0] = 1.0
    RS[0,1] = 0.0
    RS[0,2] = 0.0

    RS[1,0] = 0.0
    RS[1,1] = 1.0
    RS[1,2] = 0.0
	
    RS[2,0] = 0.0
    RS[2,1] = 0.0
    RS[2,2] = 1.0

elif (crys==1 and s==2):

    RS[0,0] =-0.5
    RS[0,1] = 0.5
    RS[0,2] = 0.5

    RS[1,0] = 0.5
    RS[1,1] =-0.5
    RS[1,2] = 0.5

    RS[2,0] = 0.5
    RS[2,1] = 0.5
    RS[2,2] =-0.5

elif (crys==1 and s==3):

    RS[0,0] = 0.5
    RS[0,1] = 0.5
    RS[0,2] = 0.0

    RS[1,0] = 0.5
    RS[1,1] = 0.0
    RS[1,2] = 0.5

    RS[2,0] = 0.0
    RS[2,1] = 0.5
    RS[2,2] = 0.5

elif (crys==2 and s==1):
 
    RS[0,0] = 1.0
    RS[0,1] = 0.0
    RS[0,2] = 0.0

    RS[1,0] = 0.0
    RS[1,1] = 1.0
    RS[1,2] = 0.0

    RS[2,0] = 0.0
    RS[2,1] = 0.0
    RS[2,2] = ca

elif (crys==2 and s==2):

    RS[0,0] =-0.5
    RS[0,1] = 0.5
    RS[0,2] = ca/2.0

    RS[1,0] = 0.5
    RS[1,1] =-0.5
    RS[1,2] = ca/2.0

    RS[2,0] = 0.5
    RS[2,1] = 0.5
    RS[2,2] =-ca/2.0

elif (crys==3 and s==1):

    RS[0,0] = 1.0
    RS[0,1] = 0.0
    RS[0,2] = 0.0

    RS[1,0] = 0.0
    RS[1,1] = ba
    RS[1,2] = 0.0

    RS[2,0] = 0.0
    RS[2,1] = 0.0
    RS[2,2] = ca

elif (crys==3 and s==2):

    RS[0,0] =-0.5
    RS[0,1] = 0.5*ba
    RS[0,2] = 0.5*ca

    RS[1,0] = 0.5
    RS[1,1] =-0.5*ba
    RS[1,2] = 0.5*ca

    RS[2,0] = 0.5
    RS[2,1] = 0.5*ba
    RS[2,2] =-0.5*ca

elif (crys==3 and s==3):

    RS[0,0] = 0.5
    RS[0,1] = 0.5*ba
    RS[0,2] = 0.0

    RS[1,0] = 0.5
    RS[1,1] = 0.0
    RS[1,2] = 0.5*ca

    RS[2,0] = 0.0
    RS[2,1] = 0.5*ba
    RS[2,2] = 0.5*ca

elif (crys==3 and s==4):

    RS[0,0] = 0.5
    RS[0,1] = 0.5*ba
    RS[0,2] = 0.0

    RS[1,0] =-0.5
    RS[1,1] = 0.5*ba
    RS[1,2] = 0.0

    RS[2,0] = 0.0
    RS[2,1] = 0.0
    RS[2,2] = ca

elif(crys==4 and s==1):

    RS[0,0] = 1.0
    RS[0,1] = 0.0
    RS[0,2] = 0.0

    RS[1,0] = ba*np.cos(gamma)
    RS[1,1] = ba*np.sin(gamma)
    RS[1,2] = 0.0

    RS[2,0] = 0.0
    RS[2,1] = 0.0
    RS[2,2] = ca

elif(crys==4 and s==2):# NEEDS TO BE checked

    RS[0,0] = 0.5
    RS[0,1] = 0.0
    RS[0,2] = -ca/2.0

    RS[1,0] = ba*np.cos(gamma)
    RS[1,1] = ba*np.sin(gamma)
    RS[1,2] = 0.0

    RS[2,0] = 0.5
    RS[2,1] = 0.0
    RS[2,2] = ca/2.0

elif (crys==5):

    RS[0,0] = 1.0
    RS[0,1] = 0.0
    RS[0,2] = 0.0

    RS[1,0] = ba*np.cos(gamma)
    RS[1,1] = ba*np.sin(gamma)
    RS[1,2] = 0.0

    RS[2,0] = ca*np.cos(beta)
    RS[2,1] = ca*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
    RS[2,2] = ca*np.sqrt(1.0+2.0*np.cos(alpha)*np.cos(beta)*np.cos(gamma)-np.cos(alpha)**2-np.cos(beta)**2-np.cos(gamma)**2)/np.sin(gamma)

else:
    RS[0,0] = 1.0
    RS[0,1] = 0.0
    RS[0,2] = 0.0

    RS[1,0] =-0.5
    RS[1,1] = np.sqrt(3.0)/2.0
    RS[1,2] = 0.0

    RS[2,0] = 0.0
    RS[2,1] = 0.0
    RS[2,2] = ca

print
print ('r1= ',RS[0,0], RS[0,1], RS[0,2])
print ('r2= ',RS[1,0], RS[1,1], RS[1,2])
print ('r3= ',RS[2,0], RS[2,1], RS[2,2])
print

GS = np.zeros(shape=(3,3))
# Generate reciprocal spanning vectors in units of 2pi/a
V = np.abs(RS[0,0]*(RS[1,1]*RS[2,2] - RS[2,1]*RS[1,2])+RS[0,1]*(RS[2,0]*RS[1,2] - RS[1,0]*RS[2,2]) + RS[0,2]*(RS[1,0]*RS[2,1] - RS[2,0]*RS[1,1]))
print ('Unit cell volume= ',V)
GS[0,0] = (RS[0,1]*RS[1,2] - RS[1,1]*RS[0,2])/V
GS[0,1] = (RS[1,0]*RS[0,2] - RS[1,2]*RS[0,0])/V
GS[0,2] = (RS[0,0]*RS[1,1] - RS[1,0]*RS[0,1])/V

GS[1,0] = (RS[1,1]*RS[2,2] - RS[2,1]*RS[1,2])/V
GS[1,1] = (RS[2,0]*RS[1,2] - RS[1,0]*RS[2,2])/V
GS[1,2] = (RS[1,0]*RS[2,1] - RS[2,0]*RS[1,1])/V

GS[2,0] = (RS[2,1]*RS[0,2] - RS[0,1]*RS[2,2])/V
GS[2,1] = (RS[0,0]*RS[2,2] - RS[2,0]*RS[0,2])/V
GS[2,2] = (RS[2,0]*RS[0,1] - RS[0,0]*RS[2,1])/V
print ()
print ('b1= ',GS[0,0], GS[0,1], GS[0,2])
print ('b2= ',GS[1,0], GS[1,1], GS[1,2])
print ('b3= ',GS[2,0], GS[2,1], GS[2,2])
print ()

# set cutoff radius squared and generate Gvect
GC = int(input('Enter Gcut^2  '))

# define reciprocal lattice spanning vectors
# here it is FCC

t1=time.time()
# estimate number of first spanning vector from gamma to cutoff radius
# and allocate and initialize G vectors:generate ordered G vectors
GGS1 = GS[0,0]**2 + GS[0,1]**2 + GS[0,2]**2
GGS2 = GS[1,0]**2 + GS[1,1]**2 + GS[1,2]**2
GGS3 = GS[2,0]**2 + GS[2,1]**2 + GS[2,2]**2
GGS = min(GGS1,GGS2,GGS3)
nr=int(np.sqrt(GC/GGS)) + 1
G = gvect(GS,GC,nr)
ll=len(G)
print("length G ",ll)
# save G vectors to a file
#gdata = open('G_VECTORS.dat','w')
#for i in range(ll):
#    print >> gdata, np.linalg.norm(G[i])**2, G[i]
#gdata.close()

# Define GRID based on parallel-piped based on b1,b2,b3
m1 = int(input('enter how many multiples of G1 '))
m2 = int(input('enter how many multiples of G2 '))
m3 = int(input('enter how many multiples of G3 '))
N1 = int(input('Enter number of grid points along 1 direction '))
N2 = int(input('Enter number of grid points along 2 direction '))
N3 = int(input('Enter number of grid points along 3 direction '))
BZ = int(input('Enter BZ '))
# open files to store BZ data
if (os.path.exists('BZ_ZONE.dat')):
    os.system("rm -f BZ_ZONE.dat")

# call procedure to create BZ's
ws.gen_kp(m1,m2,m3,N1,N2,N3,BZ,GS,G,ll)
print ()
t2=time.time()
print ('time: ','%E' %(t2-t1), 'sec')

out_df = pd.read_csv("./BZ_ZONE.dat",header=None,delim_whitespace=True)
os.system("rm -f BZ_ZONE.dat")
print('file read')
out_df.columns = ['x','y','z','radius','count']
print(out_df.head())

print(out_df.describe())
from mayavi import mlab

x = out_df[['x']].values
y = out_df[['y']].values
z = out_df[['z']].values
rad = out_df[['radius']].values

src = mlab.pipeline.scalar_scatter(x, y, z, rad)
pts = mlab.pipeline.glyph(src, scale_mode='none',scale_factor=.1)

mlab.show()