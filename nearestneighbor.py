import numpy as np
from numpy.ma.extras import ndenumerate
from pymatgen.core import Lattice, IStructure, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from math import sqrt, pi
import matplotlib.pyplot as plt


def nearestneighborfinder(points):
    """Find the nearest neighbors for a set of 2D points.

    This function calculates the pairwise distances between all points and identifies
    the nearest neighbor(s) for each point. Points are considered nearest neighbors
    if they are within 0.01 units of the minimum distance.

    Parameters
    ----------
    points : numpy.ndarray
        A 2D array of shape (n, 2) containing the x,y coordinates of n points.

    Returns
    -------
    numpy.ndarray
        A 2D array of shape (m, 4) containing the nearest neighbor pairs.
        Each row contains [x1, y1, x2, y2] where (x1,y1) and (x2,y2) are
        the coordinates of a point and its nearest neighbor.

    Notes
    -----
    - Self-neighbors are excluded by setting diagonal distances to 9999
    - Multiple nearest neighbors are possible if they are equidistant
    - The output may contain duplicate pairs in different orderings
    - If a point has fewer than 3 primary nearest neighbors, the function will look for secondary nearest neighbors.
      Secondary neighbors are those that are not already neighbors of the primary neighbors.
    """

    # Calculate pairwise distances between all points
    distances = np.zeros([points.shape[0], points.shape[0]])
    for idx, val in ndenumerate(distances):
        dist = np.linalg.norm(points[idx[0], :] - points[idx[1], :])
        distances[idx[0], idx[1]] = dist
    
    # Exclude self-neighbors
    np.fill_diagonal(distances, 9999)

    nearestneighbors = []
    primary_neighbors = []

    # First pass: connect only primary nearest neighbors
    for i in range(points.shape[0]):
        mindist = float(np.min(distances[i, :]))
        nearestpts = np.where(distances[i, :] - mindist <= 0.001)[0]
        thisptpn = []
        for n in nearestpts:
            thisptpn.append(n)
            nearestneighbors.append([points[i, 0], points[i, 1], 
                                    points[n, 0], points[n, 1]])
        primary_neighbors.append(thisptpn)

    # Second pass: only connect secondary neighbors if they don't share any primary neighbors
    # we need a way to ensure the 3rd bonds retain the symmetry
    for i in range(len(primary_neighbors)):
        # If we have fewer than 3 primary neighbors, look for secondary neighbors
        total_neighbors = primary_neighbors[i].copy()
        count= 0
        while len(total_neighbors) < 3 and count <= points.shape[0]-1:
            # Temporarily exclude primary neighbors
            distances[i, total_neighbors] = 9999
            newmin = np.min(distances[i, :])
 
            if newmin/mindist > 1.5:
                break

            secondary_neighbors = np.where(distances[i, :] - newmin <= 0.001)[0] 
            for n in secondary_neighbors:
                shared_primary = False
                for m in primary_neighbors[i]:
                    if m in primary_neighbors[n]:
                        shared_primary = True
                        break

                # Connect if shares no primary neighbors
                if shared_primary == False:
                    nearestneighbors.append([points[i, 0], points[i, 1],
                                                points[n, 0], points[n, 1]])
                    total_neighbors.append(n)
                distances[i,n] = 9999
        count+=1

    return np.array(nearestneighbors)#returns all of the points and their nearest neighbors

def check_cross(e1, e2, e3, e4):
    v1 = e2-e1
    v2 = e4-e3
    A = np.array([v1,v2])
    b = e3-e1
    try:
        x = np.linalg.solve(A,b)
        if x[0] > 0 and x[0] < 1 and x[1] > 0 and x[1] < 1:
            return True
        else:
            return False
    except:
        return False

if __name__ == "__main__":
    #testing
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
