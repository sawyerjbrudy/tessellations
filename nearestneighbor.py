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

    # Keep track of which points are already connected
    connected = set()
    nearestneighbors = []

    # First pass: connect only primary nearest neighbors
    for i in range(points.shape[0]):
        mindist = np.min(distances[i, :])
        nearestpts = np.where(distances[i, :] - mindist <= 0.01)[0]
        
        for n in nearestpts:
            pair = tuple(sorted([i, n]))
            if pair not in connected and not np.array_equal(points[n], points[i]):
                connected.add(pair)
                nearestneighbors.append([points[i, 0], points[i, 1], 
                                      points[n, 0], points[n, 1]])

    # Second pass: only connect secondary neighbors if they don't share any primary neighbors
    for i in range(points.shape[0]):
        # Find primary neighbors for this point
        mindist = np.min(distances[i, :])
        primary_neighbors = np.where(distances[i, :] - mindist <= 0.01)[0]
        total_neighbors = primary_neighbors.copy().tolist()
        
        count=0
        # If we have fewer than 3 primary neighbors, look for secondary neighbors
        while len(total_neighbors) < 3 and count <= points.shape[0]-1:
            # Temporarily exclude primary neighbors
            distances[i, primary_neighbors] = 9999
            newmin = np.min(distances[i, :])
            secondary_neighbors = np.where(distances[i, :] - newmin <= 0.01)[0] 

            for n in secondary_neighbors:
                pair = tuple(sorted([i, n]))
                if pair not in connected and not np.array_equal(points[n], points[i]):
                    # Check if this secondary neighbor shares any primary neighbors with our point
                    shares_primary = False
                    for m in primary_neighbors:
                        # Get primary neighbors of our primary neighbor
                        m_mindist = np.min(distances[m, :])
                        m_primary = np.where(distances[m, :] - m_mindist <= 0.01)[0]
                        if n in m_primary:
                            shares_primary = True
                            distances[i, n] = 9999

                            break
                    
                    # Also check if this secondary neighbor is already connected to any of our primary neighbors
                    already_connected = False
                    for m in primary_neighbors:
                        if tuple(sorted([m, n])) in connected:
                            already_connected = True
                            distances[i, n] = 9999
                            break
                    
                    if not shares_primary and not already_connected:
                        connected.add(pair)
                        nearestneighbors.append([points[i, 0], points[i, 1],
                                              points[n, 0], points[n, 1]])
                        total_neighbors.append(n)
                count+=1

    return np.array(nearestneighbors)#returns all of the points and their nearest neighbors

"""
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
"""