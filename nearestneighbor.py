import numpy as np
from numpy.ma.extras import ndenumerate
from pymatgen.core import Lattice, IStructure, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from math import sqrt, pi
import matplotlib.pyplot as plt

def nearestneighborfinder(points):
    """Function which finds the nearest neighbor for a set of points."""

    distances = np.zeros([points.shape[0],points.shape[0]])#preallocates a adjacency matrix type matrix of the distances between the points
    for idx,val in ndenumerate(distances):#enumerates through the distance matrix, getting pairs of points
        dist = np.linalg.norm(points[idx[0],:]-points[idx[1],:])#calculates the distance between each pair

        distances[idx[0],idx[1]] = dist#assigns the distance point with its corresponding distance
    np.fill_diagonal(distances, 9999) #puts an arbitrary big number along the diagonal so the closest neighbor of a point is not itself


    nearestneighbors = []
    for i in range(points.shape[0]):#loops through all points
        mindist = np.min(distances[i,:])#finds the minimum distance between the point in question and all others
        nearestpts = np.where(distances[i,:]-mindist <= .01)#finds the points corresponding to that minimum distance
        for n in range(len(nearestpts[0])):#loops through the closest points
            if points[nearestpts[0][n],0] != points[i,0] or points[nearestpts[0][n],1] != points[i,1]:#checks that it is not itself
                nearestneighbors.append([points[i,0],points[i,1],points[nearestpts[0][n],0],points[nearestpts[0][n],1]])#appends a 1x4 array with the ordered pair of both points
    return np.array(nearestneighbors)#returns all of the points and their nearest neighbors

"""
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

print(nearestneighborfinder(points))
"""