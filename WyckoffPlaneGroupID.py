import numpy as np
from numpy.ma.extras import ndenumerate
from pymatgen.core import Lattice, IStructure, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from math import sqrt, pi
import matplotlib.pyplot as plt
from nearestneighbor import nearestneighborfinder

def carttofrac(points,latticevectors):
    # turning cartesian coords into fractional coords through change of basis
    # i think it works for 3d but idk for sure
    #Outputs an array of fractional coordinates

    A = np.array([latticevectors[0], latticevectors[1], latticevectors[2]])#lattice vector array
    At = A.T#transposes

    #adds an extra column of zeroes to make 2d coord arrays nx3
    if points[0].size < 3:
        points = np.column_stack((points, np.zeros(points.shape[0])))
    ptT = points.T #transposes into 3xn

    fc = np.linalg.solve(At, ptT)#solves linear equation for fractional coordinates
    return fc.T#returns nx3 array of fractional coordinates

def twodstructureplotter(latticevectors, points, size, latticeshown=False, speciesshown=False, bondshown=False):
    """This function takes a set of cartesian coordinates and lattice vectors that represent the unit cell of your structure
    Also takes a size argument which gives the square dimensions of your structure, ie. size=5 -> 5 unit cellx5 unit cell supercell
    The three boolean args tell the function what you would like to see in the plot"""

    plt.figure()#initiates matplotlib plot

    latticecoords = np.zeros([(size+1)**2, 2])#preallocates matrix to be filled with lattice points
    n = 0#variable to enumerate through the rows of the lattice coords mat
    for i in range(size+1):#iterates over one vector
        for j in range(size+1):#iterates over other vector
            posvect = ((i*latticevectors[0]) + (j*latticevectors[1]))#defines lattice points as linear combo of lattice vectors
            latticecoords[n] = [posvect[0], posvect[1]]#assigns these points to the matrix
            n += 1

    speciescoords = np.zeros([len(points)*(size**2),2])#preallocates species coordinate matrix
    modpoints = np.column_stack((points,np.zeros(len(points))))#makes the 2d points 3d so it works with 3d lattice vectors
    k = 0#iterates through the structure for each point
    for pt in modpoints:
        for i in range(size):
            for j in range(size):
                newpointv = pt + ((i * latticevectors[0]) + (j * latticevectors[1]))#moves each point by a linear combination of lattice vects
                speciescoords[k] = [newpointv[0], newpointv[1]]#assigns these new positions in matrix
                k += 1

    lineends = nearestneighborfinder(speciescoords)#uses my nearest neighbor algorithm to find the nearest neighbors of every point

    if latticeshown == True:
        plt.plot(latticecoords[:, 0], latticecoords[:, 1], 'ob', zorder=1)#graphs lattice points
    if speciesshown == True:
        plt.plot(speciescoords[:,0], speciescoords[:,1], 'ro', markersize=3, zorder=3)#graphs actual structure
    if bondshown == True and (size > 1 or len(points) > 1):
       for line in lineends:
            x1,y1,x2,y2 = line

            plt.plot([x1,x2],[y1,y2], color = 'green')#shows bonds

    plt.axis('equal')
    plt.show()


points = np.array(
    [
    [0,0],
    [.5,.5]
    ]
)

#print(nearestneighborfinder(points,1))


#species type
species = ['C']*len(points)

#unit cell description
#square for now
a1 = np.array([1,0,0])
a2 = np.array([0,1,0])
lattice_vectors = np.array([a1, a2, [0, 0, 10]])

"""Using PyMatGen to convert the points and lattice vectors to PyMatGen Structure and analyze the space group and wyckoffs. 
My project is mostly in 2D so some of the PyMatGen stuff is a little incorrect, which is why I started making the structures
myself using the plane groups and occupied wyckoff sites"""
lattice = Lattice(lattice_vectors)
structure = Structure(lattice, species, carttofrac(points,lattice_vectors))

analyzer = SpacegroupAnalyzer(structure, symprec=1e-3)
spacegroup = analyzer.get_space_group_symbol()
symmetricalstructure = analyzer.get_symmetrized_structure()


###
#print("Space group:", spacegroup)
#print(structure)
#print(symmetricalstructure.wyckoff_letters)
#print(symmetricalstructure.wyckoff_symbols)

#twodstructureplotter(lattice_vectors,points,5,latticeshown=False, speciesshown=True,bondshown=True)