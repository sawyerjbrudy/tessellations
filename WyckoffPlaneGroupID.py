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

    A = np.array([latticevectors[0], latticevectors[1], latticevectors[2]])
    At = A.T

    if points[0].size < 3:
        points = np.column_stack((points, np.zeros(points.shape[0])))
    ptT = points.T

    fc = np.linalg.solve(At, ptT)
    return fc.T

def twodstructureplotter(latticevectors, points, size, latticeshown=False, speciesshown=False, bondshown=False):
    plt.figure()

    latticecoords = np.zeros([(size+1)**2, 2])
    n = 0
    for i in range(size+1):
        for j in range(size+1):
            posvect = ((i*latticevectors[0]) + (j*latticevectors[1]))
            latticecoords[n] = [posvect[0], posvect[1]]
            n += 1

    speciescoords = np.zeros([len(points)*(size**2),2])
    modpoints = np.column_stack((points,np.zeros(len(points))))
    k = 0
    for pt in modpoints:
        for i in range(size):
            for j in range(size):
                newpointv = pt + ((i * latticevectors[0]) + (j * latticevectors[1]))
                speciescoords[k] = [newpointv[0], newpointv[1]]
                k += 1

    lineends = nearestneighborfinder(speciescoords, 2)

    if latticeshown == True:
        plt.plot(latticecoords[:, 0], latticecoords[:, 1], 'ob', zorder=1)
    if speciesshown == True:
        plt.plot(speciescoords[:,0], speciescoords[:,1], 'ro', markersize=3, zorder=3)
    if bondshown == True and (size > 1 or len(points) > 1):
       for line in lineends:
            x1,y1,x2,y2 = line

            plt.plot([x1,x2],[y1,y2], color = 'green')

    plt.axis('equal')
    plt.show()


points = np.array(
    [
    [0,0],
    [.5,.5]
    ]
)

#print(nearestneighborfinder(points,1))


#species which makes up motif?
species = ['C']*len(points)

#unit cell description
#square for now
a1 = np.array([1,0,0])
a2 = np.array([0,1,0])
lattice_vectors = np.array([a1, a2, [0, 0, 10]])

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