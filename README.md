# Visualising_BZ
Plot the higher order Brillouin zones of Crystals

# python environment

Need a python3 environment with numpy, pandas and for graphics install mayavi

# gfortran

You also need gfortran to create the library to find k_points

The python script will see if the ws library has been created
Else it will run f2py to compile bzones_lib_counting.f90 to make the library

# Running the code

export OMP_NUM_THREADS=8 or whatever for openMP

python BZ_ZONES_counting.py

        The code will check if the ws library ws.cpython-37m-darwin.so exits
        If not it runs f2py on bzones_lib_counting.f90 to create the library

Enter parameters (there are defaults just hit enter):

choose lattice type

choose Gcut like 81

choose grid in terms of b1,b2,b3 depends on the bz order. 2x2x2 are reasonable values for most bz

Enter number of grid points for each direction: Depends on your memory but an upper limit of 400x400x400 ok

Choose bz to save

It will create a file of k_points in that bz then load them as a pandas DataFrame 
 
mayavi creates the visualisation