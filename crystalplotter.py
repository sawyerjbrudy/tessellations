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

def twodstructureplotter(latticevectors, fraccoords, size, latticeshown=False, speciesshown=False, bondshown=False, plane="",wyckoffs=[]):
    """This function takes a set of cartesian coordinates and lattice vectors that represent the unit cell of your structure
    Also takes a size argument which gives the square dimensions of your structure, ie. size=5 -> 5 unit cellx5 unit cell supercell
    The three boolean args tell the function what you would like to see in the plot"""
    
    points = np.zeros_like(fraccoords)#turns form fractional into cartesian coords
    for i in range(fraccoords.shape[0]):
        points[i,:] = latticevectors[0,:2]*fraccoords[i,0] + latticevectors[1,:2]*fraccoords[i,1]


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
    
    if plane != "" and wyckoffs != []:
        plt.title("2D Structure of "+plane+" "+str(wyckoffs))

    plt.axis("equal")
    plt.show()

def supercellmaker(latticevectors,fraccoords):
    #makes a 7x7 supercell in order to analyze the bonding in the middle of the unit cell

    points = np.zeros_like(fraccoords)#turns form fractional into cartesian coords
    for i in range(fraccoords.shape[0]):
        points[i,:] = latticevectors[0,:2]*fraccoords[i,0] + latticevectors[1,:2]*fraccoords[i,1]

    supercell = np.zeros([7,7,points.shape[0],2])
    for i in range(7):
        for j in range(7):
            supercell[i,j,:] = points + ((i*latticevectors[0,0:2]) + (j*latticevectors[1,0:2]))
    return supercell



points = np.array(
    [
    [.3333,.3333],
    [.6667,.6667]
    ]
)

lattice_vectors = np.array([[1,0,0],[.5,sqrt(3)/2,0],[0,0,1]])

if __name__ == "__main__":
    #print("Space group:", spacegroup)
    #print(structure)
    #print(symmetricalstructure.wyckoff_letters)
    #print(symmetricalstructure.wyckoff_symbols)

    twodstructureplotter(lattice_vectors,points,9,latticeshown=False, speciesshown=True,bondshown=True)