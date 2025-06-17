import matplotlib.pyplot as plt
import numpy as np
from WyckoffPlaneGroupID import carttofrac, twodstructureplotter
from math import pi, sqrt

def builda2dcrystal(planegroup : str, occupiedwyckoffpositions : list[str]) -> np.ndarray:
    """Generate fractional coordinates for a 2D crystal structure based on plane group and Wyckoff positions.

    This function takes a plane group and occupied Wyckoff positions (Wyckoff letters) and returns
    a set of fractional coordinates that represent the crystal structure. For symmetry groups with
    mirror planes, the function assigns a general location on the symmetry element for the point.

    Parameters
    ----------
    planegroup : str
        The plane group symbol (e.g. "p1", "p2", "pm", etc.)
    occupiedwyckoffpositions : list of str
        List of Wyckoff position letters that are occupied in the structure

    Returns
    -------
    numpy.ndarray
        Array of fractional coordinates for all equivalent positions in the unit cell

    Notes
    -----
    The function currently has limited control over the exact positions of points on symmetry elements.
    The implementation is quite tedious due to the need to handle all possible plane groups and
    Wyckoff positions.
    """

    #takes a plane group and occupied wyckoff positions (wyckoff letters) and gives a generalized sets of fractional points that represent the crystal
    #On symmetry groups such as mirror planes function assigns a general location on the symmetry element for the point, no control yet
    #very tedious

    #list of plane groups
    listofplanegroups = ["p1","p2","pm","pg","cm","pmm","pmg","pgg","cmm","p4","p4m","p4g","p3","p3m","p3m1","p6","p6m"]
    fraccoords = []#fractional coord list to be appended to
    if planegroup in listofplanegroups: #finds the wanted plane groups, appends the points of the wanted wps
        if planegroup == listofplanegroups[0]:
            for wp in occupiedwyckoffpositions:
                if wp == "a":
                    fraccoords.append([.5,.5])
        elif planegroup == listofplanegroups[1]:
            for wp in occupiedwyckoffpositions:
                if wp == "e":
                    fraccoords.append([.25,.25])
                    fraccoords.append([-.25,-.25])
                elif wp == "d":
                    fraccoords.append([.5,.5])
                elif wp == "c":
                    fraccoords.append([.5,0])
                elif wp == "b":
                    fraccoords.append([0,.5])
                elif wp == "a":
                    fraccoords.append([0,0])
        elif planegroup == listofplanegroups[2]:
            for wp in occupiedwyckoffpositions:
                if wp == "c":
                    fraccoords.append([.25,.25])
                    fraccoords.append([-.25,.25])
                elif wp == "b":
                    fraccoords.append([.5, .5])
                elif wp == "a":
                    fraccoords.append([0,.5])
        elif planegroup == listofplanegroups[3]:
            for wp in occupiedwyckoffpositions:
                if wp == "a":
                    fraccoords.append([.25,.25])
                    fraccoords.append([-.25,.75])
        elif planegroup == listofplanegroups[4]:
            for wp in occupiedwyckoffpositions:
                if wp == "b":
                    fraccoords.append([.25,.25])
                    fraccoords.append([-.25,.25])
                elif wp == "a":
                    fraccoords.append([0,.5])
        elif planegroup == listofplanegroups[5]:
            for wp in occupiedwyckoffpositions:
                if wp == "i":
                    fraccoords.append([.25,.25])
                    fraccoords.append([-.25,-.25])
                    fraccoords.append([-.25,.25])
                    fraccoords.append([.25,-.25])
                elif wp == "h":
                    fraccoords.append([.5,.75])
                    fraccoords.append([.5,-.75])
                elif wp == "g":
                    fraccoords.append([0,.75])
                    fraccoords.append([0,-.75])
                elif wp == "f":
                    fraccoords.append([.75,.5])
                    fraccoords.append([-.75,.5])
                elif wp == "e":
                    fraccoords.append([.75,0])
                    fraccoords.append([-.75,0])
                elif wp == "d":
                    fraccoords.append([.5,.5])
                elif wp == "c":
                    fraccoords.append([.5,0])
                elif wp == "b":
                    fraccoords.append([0,.5])
                elif wp == "a":
                    fraccoords.append([0,0])
        elif planegroup == listofplanegroups[6]:
            for wp in occupiedwyckoffpositions:
                if wp == "d":
                    fraccoords.append([.125,.125])
                    fraccoords.append([.375,.125])
                    fraccoords.append([-.125,-.125])
                    fraccoords.append([.625,-.125])
                elif wp == "c":
                    fraccoords.append([.25,.25])
                    fraccoords.append([.75,-.25])
                elif wp == "b":
                    fraccoords.append([0,.5])
                    fraccoords.append([.5,.5])
                elif wp == "a":
                    fraccoords.append([0,0])
                    fraccoords.append([.5,0])
        elif planegroup == listofplanegroups[7]:
            for wp in occupiedwyckoffpositions:
                if wp == "c":
                    fraccoords.append([.125,.125])
                    fraccoords.append([-.125,-.125])
                    fraccoords.append([.375,.625])
                    fraccoords.append([.625,.375])
                elif wp == "b":
                    fraccoords.append([0,.5])
                    fraccoords.append([.5,0])
                elif wp == "a":
                    fraccoords.append([0,0])
                    fraccoords.append([.5,.5])
        elif planegroup == listofplanegroups[8]:
            for wp in occupiedwyckoffpositions:
                if wp == "f":
                    fraccoords.append([.125,.125])
                    fraccoords.append([-.125,-.125])
                    fraccoords.append([-.125,.125])
                    fraccoords.append([.125,-.125])
                    fraccoords.append([.625,.625])
                    fraccoords.append([-.625,-.625])
                    fraccoords.append([.625,-.625])
                    fraccoords.append([-.625,.625])
                elif wp == "e":
                    fraccoords.append([0,.29])
                    fraccoords.append([0,-.29])
                    fraccoords.append([0, .79])
                    fraccoords.append([0, -.79])
                elif wp == "d":
                    fraccoords.append([.25,0])
                    fraccoords.append([-.25,0])
                    fraccoords.append([-.75,0])
                    fraccoords.append([.75,0])
                elif wp == "c":
                    fraccoords.append([.29,.29])
                    fraccoords.append([.79,.29])
                    fraccoords.append([.79, .79])
                    fraccoords.append([1.29, .79])
                elif wp == "b":
                    fraccoords.append([0,.5])
                    fraccoords.append([.5,1])
                elif wp == "a":
                    fraccoords.append([0,0])
                    fraccoords.append([.5,.5])
        elif planegroup == listofplanegroups[9]:
            for wp in occupiedwyckoffpositions:
                if wp == "d":
                    fraccoords.append([.125,.125])
                    fraccoords.append([-.125,-.125])
                    fraccoords.append([-.125,.125])
                    fraccoords.append([.125,-.125])
                    fraccoords.append([.625,.625])
                    fraccoords.append([-.625,-.625])
                    fraccoords.append([.625,-.625])
                    fraccoords.append([-.625,.625])
                elif wp == "c":
                    fraccoords.append([.5,0])
                    fraccoords.append([0,.5])
                elif wp == "b":
                    fraccoords.append([.5,.5])
                elif wp == "a":
                    fraccoords.append([0,0])
        elif planegroup == listofplanegroups[10]:
            for wp in occupiedwyckoffpositions:
                if wp == "g":
                    fraccoords.append([.25,.75])
                    fraccoords.append([-.25,.75])
                    fraccoords.append([-.25,-.75])
                    fraccoords.append([.25,-.75])
                    fraccoords.append([.75,.25])
                    fraccoords.append([-.75,.25])
                    fraccoords.append([.75,-.25])
                    fraccoords.append([-.75,-.25])
                elif wp == "f":
                    fraccoords.append([.25,.25])
                    fraccoords.append([-.25,.25])
                    fraccoords.append([-.25,-.25])
                    fraccoords.append([.25,-.25])
                elif wp == "e":
                    fraccoords.append([.25,.5])
                    fraccoords.append([-.25,.5])
                    fraccoords.append([.5,.25])
                    fraccoords.append([.5,-.25])
                elif wp == "d":
                    fraccoords.append([.25,0])
                    fraccoords.append([-.25,0])
                    fraccoords.append([0,-.25])
                    fraccoords.append([0, .25])
                elif wp == "c":
                    fraccoords.append([.5,0])
                    fraccoords.append([0,.5])
                elif wp == "b":
                    fraccoords.append([.5,.5])
                elif wp == "a":
                    fraccoords.append([0,0])
        elif planegroup == listofplanegroups[11]:
            for wp in occupiedwyckoffpositions:
                if wp == "d":
                    fraccoords.append([.125,.75])
                    fraccoords.append([.75,-.125])
                    fraccoords.append([-.125,-.75])
                    fraccoords.append([-.75,.125])
                    fraccoords.append([.375,1.25])
                    fraccoords.append([.625,-.25])
                    fraccoords.append([1.25,.625])
                    fraccoords.append([-.25,-.375])
                elif wp == "c":
                    fraccoords.append([.25,.75])
                    fraccoords.append([-.25,.25])
                    fraccoords.append([.25,.25])
                    fraccoords.append([.75,-.25])
                elif wp == "b":
                    fraccoords.append([0,.5])
                    fraccoords.append([.5,0])
                elif wp == "a":
                    fraccoords.append([0,0])
                    fraccoords.append([.5,.5])
        elif planegroup == listofplanegroups[12]:
            for wp in occupiedwyckoffpositions:
                if wp == "d":
                    fraccoords.append([.125,.75])
                    fraccoords.append([-.75,-.625])
                    fraccoords.append([.625,-.125])
                elif wp == "c":
                    fraccoords.append([.6667,.3333])
                elif wp == "b":
                    fraccoords.append([.3333, .6667])
                elif wp == "a":
                    fraccoords.append([0,0])
        elif planegroup == listofplanegroups[13]:
            for wp in occupiedwyckoffpositions:
                if wp == "e":
                    fraccoords.append([.125,.75])
                    fraccoords.append([-.75,-.625])
                    fraccoords.append([.625,-.125])
                    fraccoords.append([-.75,-.125])
                    fraccoords.append([.625,.75])
                    fraccoords.append([.125, -.625])
                elif wp == "d":
                    fraccoords.append([.25,-.25])
                    fraccoords.append([.25,.5])
                    fraccoords.append([-.5,-.25])
                elif wp == "c":
                    fraccoords.append([.6667,.3333])
                elif wp == "b":
                    fraccoords.append([.3333, .6667])
                elif wp == "a":
                    fraccoords.append([0,0])
        elif planegroup == listofplanegroups[14]:
            for wp in occupiedwyckoffpositions:
                if wp == "d":
                    fraccoords.append([.125,.75])
                    fraccoords.append([-.75,-.625])
                    fraccoords.append([.625,-.125])
                    fraccoords.append([.75,.125])
                    fraccoords.append([-.625,-.75])
                    fraccoords.append([-.125,.625])
                elif wp == "c":
                    fraccoords.append([.25,0])
                    fraccoords.append([0,.25])
                    fraccoords.append([-.25,-.25])
                elif wp == "b":
                    fraccoords.append([.3333,.6667])
                    fraccoords.append([.6667,.3333])
                elif wp == "a":
                    fraccoords.append([0,0])
        elif planegroup == listofplanegroups[15]:
            for wp in occupiedwyckoffpositions:
                if wp == "d":
                    fraccoords.append([.125,.75])
                    fraccoords.append([-.75,-.625])
                    fraccoords.append([.625,-.125])
                    fraccoords.append([-.125,-.75])
                    fraccoords.append([.75,.625])
                    fraccoords.append([-.625,.125])
                elif wp == "c":
                    fraccoords.append([.5,0])
                    fraccoords.append([0,.5])
                    fraccoords.append([.5,.5])
                elif wp == "b":
                    fraccoords.append([.3333,.6667])
                    fraccoords.append([.6667,.3333])
                elif wp == "a":
                    fraccoords.append([0,0])
        elif planegroup == listofplanegroups[16]:
            for wp in occupiedwyckoffpositions:
                if wp == "f":
                    fraccoords.append([.125,.75])
                    fraccoords.append([-.75,-.625])
                    fraccoords.append([.625,-.125])
                    fraccoords.append([-.125,-.75])
                    fraccoords.append([.75,.625])
                    fraccoords.append([-.625,.125])
                    fraccoords.append([-.75,-.125])
                    fraccoords.append([.625,.75])
                    fraccoords.append([.125,-.625])
                    fraccoords.append([.75,.125])
                    fraccoords.append([-.625,-.75])
                    fraccoords.append([-.125,.625])
                elif wp == "e":
                    fraccoords.append([.25,-.25])
                    fraccoords.append([.25,.5])
                    fraccoords.append([-.5,-.25])
                    fraccoords.append([-.25,.25])
                    fraccoords.append([-.25,-.5])
                    fraccoords.append([.5,.25])
                elif wp == "d":
                    fraccoords.append([.25,0])
                    fraccoords.append([0,.25])
                    fraccoords.append([-.25,-.25])
                    fraccoords.append([-.25,0])
                    fraccoords.append([0,-.25])
                    fraccoords.append([.25,.25])
                elif wp == "c":
                    fraccoords.append([.5,0])
                    fraccoords.append([0,.5])
                    fraccoords.append([.5,.5])
                elif wp == "b":
                    fraccoords.append([.3333,.6667])
                    fraccoords.append([.6667,.3333])
                elif wp == "a":
                    fraccoords.append([0,0])

    else:
        print("Planegroup not found")

    return np.array(fraccoords)

def generallatticevectors(planegroup: str):
    """Return lattice vectors for a given 2D plane group.

    Parameters
    ----------
    planegroup : str
        The plane group symbol (e.g. "p1", "p2", "pm", etc.)

    Returns
    -------
    numpy.ndarray
        A 3x3 array containing the lattice vectors. The first two rows contain
        the 2D lattice vectors, while the third row is [0,0,10] for compatibility
        with 3D operations.

    Notes
    -----
    The function maps plane groups to one of five standard lattice vector sets:
    1. p1: General oblique lattice
    2. p2, pm, pg, pgg, pmm, pmg: Rectangular lattice
    3. cm, cmm: Centered rectangular lattice
    4. p4, p4m, p4g: Square lattice
    5. p3, p3m, p3m1, p6, p6m: Hexagonal lattice
    """
    listofplanegroups = ["p1", "p2", "pm", "pg", "cm", "pmm", "pmg", "pgg", "cmm", 
                        "p4", "p4m", "p4g", "p3", "p3m", "p3m1", "p6", "p6m"]
    
    latticevectorsset = [
        np.array([[.93969,.34202,0], [.98481,.17365,0], [0,0,10]]),  # p1
        np.array([[1.5,0,0], [0,1,0], [0,0,10]]),                     # p2,pm,pg,pgg,pmm,pmg
        np.array([[.75,.5,0], [.75,-.5,0], [0,0,10]]),               # cm,cmm
        np.array([[1,0,0], [0,1,0], [0,0,10]]),                       # p4,p4m,p4g
        np.array([[.5,sqrt(3)/2,0], [1,0,0], [0,0,10]])              # p3,p3m,p3m1,p6,p6m
    ]

    if planegroup not in listofplanegroups:
        print("Planegroup not found")
        return None

    if planegroup == "p1":
        return latticevectorsset[0]
    elif planegroup in ["p2", "pm", "pg", "pgg", "pmm", "pmg"]:
        return latticevectorsset[1]
    elif planegroup in ["cm", "cmm"]:
        return latticevectorsset[2]
    elif planegroup in ["p4", "p4m", "p4g"]:
        return latticevectorsset[3]
    elif planegroup in ["p3", "p3m", "p3m1", "p6", "p6m"]:
        return latticevectorsset[4]

if __name__ == "__main__":
    #all the testing
    pts = builda2dcrystal("p6", ["a","c","d"])#fractional coords
    latticevectors = generallatticevectors("p6")#lattice vectors

    cartpts = np.zeros_like(pts)#turns form fractional into cartesian coords
    for i in range(pts.shape[0]):
        cartpts[i,:] = latticevectors[0,:2]*pts[i,0] + latticevectors[1,:2]*pts[i,1]
    print(cartpts)
    twodstructureplotter(latticevectors,cartpts,8,latticeshown=False,speciesshown=False,bondshown=True)#plots the structure