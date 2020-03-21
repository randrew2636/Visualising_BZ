# Richard Andrew 2020
# 
import numpy as np
import pandas as pd

import time
import os
#--------------------------------------------------------

# Run f2py compiling the openMP fortran file into a python-callable library

os.system("f2py3.7 -c -m ws bzones_lib_counting.f90 --f90flags=-fopenmp -lgomp")

# Import this library

import ws

# NOTE 'ws.*.so' is the library and 'ws' is the import name

#----------------------------------------------
# Gvect generator function

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
crys = int(input('Enter (default 1): ') or '1')

# Enter crystal type and associated parameters

ca=1.0
ba=1.0
if(crys==1):
    print ('P-1 I-2 F-3')
    s = int(input('Enter: (default 3) ') or '3')
elif(crys==2):
    print ('P-1 I-2')
    s = int(input('Enter: ') or '2')
    ca = float(input('enter c/a (default 0.5) ') or '0.5')
elif(crys==3):
    print ('P-1 I-2 F-3 C-4')
    s = int(input('Enter: (default 3) ') or '3')
    ba = float(input('enter b/a (default 0.5) ') or '0.5')
    ca = float(input('enter c/a (default 0.5) ') or 'o.5')
elif(crys==4):
    print ('P-1 I-2-NOT CHECKED')
    s = int(input('Enter: (default 1) ') or '1')
    ba = float(input('enter b/a (default 0.5) ') or '0.5')
    ca = float(input('enter c/a (default 0.5) ') or '0.5')
    gamma = float(input('enter gamma angle between a anc b vectors '))
    gamma=gamma*PI/180.0
elif(crys==5):
    ba = float(input('enter b/a (default 0.5) ') or '0.5')
    ca = float(input('enter c/a (default 0.5) ') or '0.5')

    alpha = float(input('enter alpha angle between b anc c vectors (default 30) ') or '30')
    beta  = float(input('enter beta  angle between a anc c vectors (default 30) ') or '30')
    gamma = float(input('enter gamma angle between a anc b vectors (default 30) ') or '30')

    alpha=alpha*PI/180.0
    beta=beta*PI/180.0
    gamma=gamma*PI/180.0
else:
    ca = float(input('enter c/a (default 0.5) ') or '0.5')

print ()

# RS are the spanning vectors

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

# Generate reciprocal spanning vectors in units of 2pi/a

GS = np.zeros(shape=(3,3))

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

GC = int(input('Enter Gcut^2 (default 25)  ') or '25')

# define reciprocal lattice spanning vectors

t1=time.time()

# estimate number of first spanning vector from gamma to cutoff radius
# and allocate and initialize G vectors:generate ordered G vectors

GGS1 = GS[0,0]**2 + GS[0,1]**2 + GS[0,2]**2
GGS2 = GS[1,0]**2 + GS[1,1]**2 + GS[1,2]**2
GGS3 = GS[2,0]**2 + GS[2,1]**2 + GS[2,2]**2
GGS = min(GGS1,GGS2,GGS3)

nr=int(np.sqrt(GC/GGS)) + 1

G = gvect(GS,GC,nr) # generate G vector with length squared less than Gcut^2

ll=len(G)
print("length G ",ll)

# Define GRID based on parallel-piped based on b1,b2,b3

m1 = int(input('enter how many multiples of G1 (default 2) ') or '2')
m2 = int(input('enter how many multiples of G2 (default 2) ') or '2')
m3 = int(input('enter how many multiples of G3 (default 2) ') or '2')
N1 = int(input('Enter number of grid points along 1 direction (default 300) ') or '300')
N2 = int(input('Enter number of grid points along 2 direction (default 300) ') or '300')
N3 = int(input('Enter number of grid points along 3 direction (default 300) ') or '300')

# choose bz

BZ = int(input('Enter BZ (default 16) ') or '16') 

# open files to store BZ data

if (os.path.exists('BZ_ZONE.dat')):
    os.system("rm -f BZ_ZONE.dat")

# call fortran procedure through ws library to create k-points within bz
# Results saved to a datafile

ws.gen_kp(m1,m2,m3,N1,N2,N3,BZ,GS,G,ll)

print ()
t2=time.time()
print ('time: ','%E' %(t2-t1), 'sec')

# read k-points into a DataFrame

out_df = pd.read_csv("./BZ_ZONE.dat",header=None,delim_whitespace=True)
os.system("rm -f BZ_ZONE.dat") # delete k-point datafile
print('file read')
out_df.columns = ['x','y','z','radius','BZ']
print(out_df.head())

print(out_df.describe())

# import visualisation library

from mayavi import mlab

mlab.figure(figure=None, bgcolor=(0.7,0.7,0.7), fgcolor=None, engine=None, size=(4000, 3500))

x = out_df[['x']].values
y = out_df[['y']].values
z = out_df[['z']].values
rad = out_df[['radius']].values

# create scatter pipeline with x,y,z and radius to define depth colors

src = mlab.pipeline.scalar_scatter(x, y, z, rad)

# plot points in a mayavi sceen window

#pts = mlab.pipeline.glyph(src, scale_mode='none',scale_factor=.1,colormap="winter")
pts = mlab.pipeline.glyph(src, scale_mode='none',scale_factor=.1,colormap="gray")

mlab.show()