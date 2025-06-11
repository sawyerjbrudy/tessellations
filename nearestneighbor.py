import numpy as np
from numpy.ma.extras import ndenumerate
from pymatgen.core import Lattice, IStructure, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from math import sqrt, pi
import matplotlib.pyplot as plt

def nearestneighborfinder(points,closeness):
    distances = np.zeros([points.shape[0],points.shape[0]])
    for idx,val in ndenumerate(distances):
        dist = np.linalg.norm(points[idx[0],:]-points[idx[1],:])

        distances[idx[0],idx[1]] = dist
    np.fill_diagonal(distances, 9999)


    nearestneighbors = []
    for n in range(closeness):
        for i in range(points.shape[0]):
            mindist = np.min(distances[i,:])
            nearestpts = np.where(distances[i,:]-mindist <= .01)
            for n in range(len(nearestpts[0])):
                if points[nearestpts[0][n],0] != points[i,0] or points[nearestpts[0][n],1] != points[i,1]:
                    nearestneighbors.append([points[i,0],points[i,1],points[nearestpts[0][n],0],points[nearestpts[0][n],1]])
    return np.array(nearestneighbors)


points = np.array(
    [
    [0,0],
    [.5,.5],
    [sqrt(2)/2,0],
    [1,1],
    [1.5,1.5],
    [sqrt(2)/2+1,1]
    ]
)

print(nearestneighborfinder(points,1))