import matplotlib.pyplot as plt
import numpy as np
from crystalplotter import carttofrac, twodstructureplotter
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
    
    occupiedwyckoffpositions = sorted(occupiedwyckoffpositions)

    #list of plane groups
    listofplanegroups = ["p1","p2","pm","pg","cm","pmm","pmg","pgg","cmm","p4","p4m","p4g","p3","p3m","p3m1","p6","p6m"]
    fraccoords = []  # fractional coord list to be appended to
    if planegroup in listofplanegroups:  # finds the wanted plane groups, appends the points of the wanted wps
        if planegroup == listofplanegroups[0]:
            for wp in occupiedwyckoffpositions:
                if wp == "a":
                    x_a = input("x_a: ")
                    y_a = input("y_a: ")

                    fraccoords.append([float(x_a), float(y_a)])
        elif planegroup == listofplanegroups[1]:
            for wp in occupiedwyckoffpositions:
                if wp == "e":
                    x_e = input("x_e: ")
                    y_e = input("y_e: ")

                    fraccoords.append([float(x_e), float(y_e)])
                    fraccoords.append([-float(x_e), -float(y_e)])
                elif wp == "d":
                    fraccoords.append([.5, .5])
                elif wp == "c":
                    fraccoords.append([.5, 0])
                elif wp == "b":
                    fraccoords.append([0, .5])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == listofplanegroups[2]:
            for wp in occupiedwyckoffpositions:
                if wp == "c":
                    x_c = input("x_c: ")
                    y_c = input("y_c: ")

                    fraccoords.append([float(x_c), float(y_c)])
                    fraccoords.append([-float(x_c), float(y_c)])
                elif wp == "b":
                    y_b = input("y_b: ")

                    fraccoords.append([.5, float(y_b)])
                elif wp == "a":
                    y_a = input("y_a: ")

                    fraccoords.append([0, float(y_a)])
        elif planegroup == listofplanegroups[3]:
            for wp in occupiedwyckoffpositions:
                if wp == "a":
                    x_a = input("x_a: ")
                    y_a = input("y_a: ")

                    fraccoords.append([float(x_a), float(y_a)])
                    fraccoords.append([-float(x_a), float(y_a) + .5])
        elif planegroup == listofplanegroups[4]:
            for wp in occupiedwyckoffpositions:
                if wp == "b":
                    x_b = input("x_b: ")
                    y_b = input("y_b: ")

                    fraccoords.append([float(x_b), float(y_b)])
                    fraccoords.append([-float(x_b), float(y_b)])
                    fraccoords.append([float(x_b) + .5, float(y_b) + .5])
                    fraccoords.append([-float(x_b) + .5, float(y_b) + .5])
                elif wp == "a":
                    y_a = input("y_a: ")

                    fraccoords.append([0, float(y_a)])
                    fraccoords.append([.5, float(y_a) + .5])
        elif planegroup == listofplanegroups[5]:
            for wp in occupiedwyckoffpositions:
                if wp == "i":
                    x_i = input("x_i: ")
                    y_i = input("y_i: ")

                    fraccoords.append([float(x_i), float(y_i)])
                    fraccoords.append([-float(x_i), -float(y_i)])
                    fraccoords.append([-float(x_i), float(y_i)])
                    fraccoords.append([float(x_i), -float(y_i)])
                elif wp == "h":
                    y_h = input("y_h: ")

                    fraccoords.append([.5, float(y_h)])
                    fraccoords.append([.5, -float(y_h)])
                elif wp == "g":
                    y_g = input("y_g: ")

                    fraccoords.append([0, float(y_g)])
                    fraccoords.append([0, -float(y_g)])
                elif wp == "f":
                    x_f = input("x_f: ")

                    fraccoords.append([float(x_f), .5])
                    fraccoords.append([-float(x_f), .5])
                elif wp == "e":
                    x_e = input("x_e: ")

                    fraccoords.append([float(x_e), 0])
                    fraccoords.append([-float(x_e), 0])
                elif wp == "d":
                    fraccoords.append([.5, .5])
                elif wp == "c":
                    fraccoords.append([.5, 0])
                elif wp == "b":
                    fraccoords.append([0, .5])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == listofplanegroups[6]:
            for wp in occupiedwyckoffpositions:
                if wp == "d":
                    x_d = input("x_d: ")
                    y_d = input("y_d: ")

                    fraccoords.append([float(x_d), float(y_d)])
                    fraccoords.append([-float(x_d), -float(y_d)])
                    fraccoords.append([-float(x_d) + .5, float(y_d)])
                    fraccoords.append([float(x_d) + .5, -float(y_d)])
                elif wp == "c":
                    y_c = input("y_c: ")

                    fraccoords.append([.25, float(y_c)])
                    fraccoords.append([.75, -float(y_c)])
                elif wp == "b":
                    fraccoords.append([0, .5])
                    fraccoords.append([.5, .5])
                elif wp == "a":
                    fraccoords.append([0, 0])
                    fraccoords.append([.5, 0])
        elif planegroup == listofplanegroups[7]:
            for wp in occupiedwyckoffpositions:
                if wp == "c":
                    x_c = input("x_c: ")
                    y_c = input("y_c: ")

                    fraccoords.append([float(x_c), float(y_c)])
                    fraccoords.append([-float(x_c), -float(y_c)])
                    fraccoords.append([-float(x_c) + .5, float(y_c) + .5])
                    fraccoords.append([float(x_c) + .5, -float(y_c) + .5])
                elif wp == "b":
                    fraccoords.append([0, .5])
                    fraccoords.append([.5, 0])
                elif wp == "a":
                    fraccoords.append([0, 0])
                    fraccoords.append([.5, .5])
        elif planegroup == listofplanegroups[8]:
            for wp in occupiedwyckoffpositions:
                if wp == "f":
                    x_f = input("x_f: ")
                    y_f = input("y_f: ")

                    fraccoords.append([float(x_f), float(y_f)])
                    fraccoords.append([-float(x_f), -float(y_f)])
                    fraccoords.append([-float(x_f), float(y_f)])
                    fraccoords.append([float(x_f), -float(y_f)])
                    fraccoords.append([float(x_f) + .5, float(y_f) + .5])
                    fraccoords.append([-float(x_f) + .5, -float(y_f) + .5])
                    fraccoords.append([-float(x_f) + .5, float(y_f) + .5])
                    fraccoords.append([float(x_f) + .5, -float(y_f) + .5])
                elif wp == "e":
                    y_e = input("y_e: ")

                    fraccoords.append([0, float(y_e)])
                    fraccoords.append([0, -float(y_e)])
                    fraccoords.append([.5, float(y_e) + .5])
                    fraccoords.append([.5, -float(y_e) + .5])
                elif wp == "d":
                    x_d = input("x_d: ")

                    fraccoords.append([float(x_d), 0])
                    fraccoords.append([-float(x_d), 0])
                    fraccoords.append([float(x_d) + .5, .5])
                    fraccoords.append([-float(x_d) + .5, .5])
                elif wp == "c":
                    fraccoords.append([.25, .25])
                    fraccoords.append([.75, .25])
                    fraccoords.append([.75, .75])
                    fraccoords.append([1.25, .75])
                elif wp == "b":
                    fraccoords.append([0, .5])
                    fraccoords.append([.5, 1])
                elif wp == "a":
                    fraccoords.append([0, 0])
                    fraccoords.append([.5, .5])
        elif planegroup == listofplanegroups[9]:
            for wp in occupiedwyckoffpositions:
                if wp == "d":
                    x_d = input("x_d: ")
                    y_d = input("y_d: ")

                    fraccoords.append([float(x_d), float(y_d)])
                    fraccoords.append([-float(x_d), -float(y_d)])
                    fraccoords.append([-float(y_d), float(x_d)])
                    fraccoords.append([float(y_d), -float(x_d)])
                elif wp == "c":
                    fraccoords.append([.5, 0])
                    fraccoords.append([0, .5])
                elif wp == "b":
                    fraccoords.append([.5, .5])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == listofplanegroups[10]:
            for wp in occupiedwyckoffpositions:
                if wp == "g":
                    x_g = input("x_g: ")
                    y_g = input("y_g: ")

                    fraccoords.append([float(x_g), float(y_g)])
                    fraccoords.append([-float(x_g), -float(y_g)])
                    fraccoords.append([-float(y_g), float(x_g)])
                    fraccoords.append([float(y_g), -float(x_g)])
                    fraccoords.append([-float(x_g), float(y_g)])
                    fraccoords.append([float(x_g), -float(y_g)])
                    fraccoords.append([float(y_g), float(x_g)])
                    fraccoords.append([-float(y_g), -float(x_g)])
                elif wp == "f":
                    x_f = input("x_f: ")

                    fraccoords.append([float(x_f), float(x_f)])
                    fraccoords.append([-float(x_f), -float(x_f)])
                    fraccoords.append([-float(x_f), float(x_f)])
                    fraccoords.append([float(x_f), -float(x_f)])
                elif wp == "e":
                    x_e = input("x_e: ")

                    fraccoords.append([float(x_e), .5])
                    fraccoords.append([-float(x_e), .5])
                    fraccoords.append([.5, float(x_e)])
                    fraccoords.append([.5, -float(x_e)])
                elif wp == "d":
                    x_d = input("x_d: ")

                    fraccoords.append([float(x_d), 0])
                    fraccoords.append([-float(x_d), 0])
                    fraccoords.append([0, float(x_d)])
                    fraccoords.append([0, -float(x_d)])
                elif wp == "c":
                    fraccoords.append([.5, 0])
                    fraccoords.append([0, .5])
                elif wp == "b":
                    fraccoords.append([.5, .5])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == listofplanegroups[11]:
            for wp in occupiedwyckoffpositions:
                if wp == "d":
                    x_d = input("x_d: ")
                    y_d = input("y_d: ")

                    fraccoords.append([float(x_d), float(y_d)])
                    fraccoords.append([-float(x_d), -float(y_d)])
                    fraccoords.append([-float(y_d), float(x_d)])
                    fraccoords.append([float(y_d), -float(x_d)])
                    fraccoords.append([-float(x_d) + .5, float(y_d) + .5])
                    fraccoords.append([float(x_d) + .5, -float(y_d) + .5])
                    fraccoords.append([float(y_d) + .5, float(x_d) + .5])
                    fraccoords.append([-float(y_d) + .5, -float(x_d) + .5])
                elif wp == "c":
                    x_c = input("x_c: ")

                    fraccoords.append([float(x_c), float(x_c) + .5])
                    fraccoords.append([-float(x_c), -float(x_c) + .5])
                    fraccoords.append([-float(x_c) + .5, float(x_c)])
                    fraccoords.append([float(x_c) + .5, -float(x_c)])
                elif wp == "b":
                    fraccoords.append([0, .5])
                    fraccoords.append([.5, 0])
                elif wp == "a":
                    fraccoords.append([0, 0])
                    fraccoords.append([.5, .5])
        elif planegroup == listofplanegroups[12]:
            for wp in occupiedwyckoffpositions:
                if wp == "d":
                    x_d = input("x_d: ")
                    y_d = input("y_d: ")

                    fraccoords.append([float(x_d), float(y_d)])
                    fraccoords.append([-float(y_d), float(x_d) - float(y_d)])
                    fraccoords.append([-float(x_d) + float(y_d), -float(x_d)])
                elif wp == "c":
                    fraccoords.append([2/3, 1/3])
                elif wp == "b":
                    fraccoords.append([1/3, 2/3])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == listofplanegroups[13]:
            for wp in occupiedwyckoffpositions:
                if wp == "e":
                    x_e = input("x_e: ")
                    y_e = input("y_e: ")

                    fraccoords.append([float(x_e), float(y_e)])
                    fraccoords.append([-float(y_e), float(x_e) - float(y_e)])
                    fraccoords.append([-float(x_e) + float(y_e), -float(x_e)])
                    fraccoords.append([-float(y_e), -float(x_e)])
                    fraccoords.append([-float(x_e) + float(y_e), float(y_e)])
                    fraccoords.append([float(x_e), float(x_e) - float(y_e)])
                elif wp == "d":
                    x_d = input("x_d: ")

                    fraccoords.append([float(x_d), -float(x_d)])
                    fraccoords.append([float(x_d), 2 * float(x_d)])
                    fraccoords.append([-2 * float(x_d), -float(x_d)])
                elif wp == "c":
                    fraccoords.append([2/3, 1/3])
                elif wp == "b":
                    fraccoords.append([1/3, 2/3])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == listofplanegroups[14]:
            for wp in occupiedwyckoffpositions:
                if wp == "d":
                    x_d = input("x_d: ")
                    y_d = input("y_d: ")

                    fraccoords.append([float(x_d), float(y_d)])
                    fraccoords.append([-float(y_d), float(x_d) - float(y_d)])
                    fraccoords.append([-float(x_d) + float(y_d), -float(x_d)])
                    fraccoords.append([float(y_d), float(x_d)])
                    fraccoords.append([float(x_d) - float(y_d), -float(y_d)])
                    fraccoords.append([-float(x_d), -float(x_d) + float(y_d)])
                elif wp == "c":
                    x_c = input("x_c: ")

                    fraccoords.append([float(x_c), 0])
                    fraccoords.append([0, float(x_c)])
                    fraccoords.append([-float(x_c), -float(x_c)])
                elif wp == "b":
                    fraccoords.append([1/3, 2/3])
                    fraccoords.append([2/3, 1/3])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == listofplanegroups[15]:
            for wp in occupiedwyckoffpositions:
                if wp == "d":
                    x_d = input("x_d: ")
                    y_d = input("y_d: ")

                    fraccoords.append([float(x_d), float(y_d)])
                    fraccoords.append([-float(y_d), float(x_d) - float(y_d)])
                    fraccoords.append([-float(x_d) + float(y_d), -float(x_d)])
                    fraccoords.append([-float(x_d), -float(y_d)])
                    fraccoords.append([float(y_d), -float(x_d) + float(y_d)])
                    fraccoords.append([float(x_d) - float(y_d), float(x_d)])
                elif wp == "c":
                    fraccoords.append([.5, 0])
                    fraccoords.append([0, .5])
                    fraccoords.append([.5, .5])
                elif wp == "b":
                    fraccoords.append([1/3, 2/3])
                    fraccoords.append([2/3, 1/3])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == listofplanegroups[16]:
            for wp in occupiedwyckoffpositions:
                if wp == "f":
                    x_f = input("x_f: ")
                    y_f = input("y_f: ")

                    fraccoords.append([float(x_f), float(y_f)])
                    fraccoords.append([-float(y_f), float(x_f) - float(y_f)])
                    fraccoords.append([-float(x_f) + float(y_f), -float(x_f)])
                    fraccoords.append([-float(x_f), -float(y_f)])
                    fraccoords.append([float(y_f), -float(x_f) + float(y_f)])
                    fraccoords.append([float(x_f) - float(y_f), float(x_f)])
                    fraccoords.append([-float(y_f), -float(x_f)])
                    fraccoords.append([-float(x_f) + float(y_f), float(y_f)])
                    fraccoords.append([float(x_f), float(x_f) - float(y_f)])
                    fraccoords.append([float(y_f), float(x_f)])
                    fraccoords.append([float(x_f) - float(y_f), -float(y_f)])
                    fraccoords.append([-float(x_f), -float(x_f) + float(y_f)])
                elif wp == "e":
                    x_e = input("x_e: ")

                    fraccoords.append([float(x_e), -float(x_e)])
                    fraccoords.append([float(x_e), 2 * float(x_e)])
                    fraccoords.append([-2 * float(x_e), -float(x_e)])
                    fraccoords.append([-float(x_e), float(x_e)])
                    fraccoords.append([-float(x_e), -2 * float(x_e)])
                    fraccoords.append([2 * float(x_e), float(x_e)])
                elif wp == "d":
                    x_d = input("x_d: ")

                    fraccoords.append([float(x_d), 0])
                    fraccoords.append([0, float(x_d)])
                    fraccoords.append([-float(x_d), -float(x_d)])
                    fraccoords.append([-float(x_d), 0])
                    fraccoords.append([0, -float(x_d)])
                    fraccoords.append([float(x_d), float(x_d)])
                elif wp == "c":
                    fraccoords.append([.5, 0])
                    fraccoords.append([0, .5])
                    fraccoords.append([.5, .5])
                elif wp == "b":
                    fraccoords.append([1/3, 2/3])
                    fraccoords.append([2/3, 1/3])
                elif wp == "a":
                    fraccoords.append([0, 0])
                elif wp == "honeycomb":
                    fraccoords.append([1/3,1/3])
                    fraccoords.append([2/3,2/3])

        return np.array(fraccoords)
    print("Planegroup not found")
    return np.array([])

def builda2dcrystal_batch(planegroup: str, wyckoff_dict: dict) -> np.ndarray:
    """Generate fractional coordinates for a 2D crystal structure based on plane group and Wyckoff positions.

    Parameters
    ----------
    planegroup : str
        The plane group symbol (e.g. "p1", "p2", "pm", etc.)
    wyckoff_dict : dict
        Dictionary where keys are Wyckoff letters and values are either a list of coordinates (for positions that need input) or None (for positions that do not)

    Returns
    -------
    numpy.ndarray
        Array of fractional coordinates for all equivalent positions in the unit cell
    """
    listofplanegroups = ["p1","p2","pm","pg","cm","pmm","pmg","pgg","cmm","p4","p4m","p4g","p3","p3m","p3m1","p6","p6m"]
    fraccoords = []
    if planegroup in listofplanegroups:
        if planegroup == "p1":
            for wp, coords in wyckoff_dict.items():
                if wp == "a" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
        elif planegroup == "p2":
            for wp, coords in wyckoff_dict.items():
                if wp == "e" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[0]), -float(coords[1])])
                elif wp == "d":
                    fraccoords.append([.5, .5])
                elif wp == "c":
                    fraccoords.append([.5, 0])
                elif wp == "b":
                    fraccoords.append([0, .5])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == "pm":
            for wp, coords in wyckoff_dict.items():
                if wp == "c" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[0]), float(coords[1])])
                elif wp == "b" and coords is not None:
                    fraccoords.append([.5, float(coords[0])])
                elif wp == "a" and coords is not None:
                    fraccoords.append([0, float(coords[0])])
        elif planegroup == "pg":
            for wp, coords in wyckoff_dict.items():
                if wp == "a" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[0]), float(coords[1]) + .5])
        elif planegroup == "cm":
            for wp, coords in wyckoff_dict.items():
                if wp == "b" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[0]), float(coords[1])])
                    fraccoords.append([float(coords[0]) + .5, float(coords[1]) + .5])
                    fraccoords.append([-float(coords[0]) + .5, float(coords[1]) + .5])
                elif wp == "a" and coords is not None:
                    fraccoords.append([0, float(coords[0])])
                    fraccoords.append([.5, float(coords[0]) + .5])
        elif planegroup == "pmm":
            for wp, coords in wyckoff_dict.items():
                if wp == "i" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[0]), -float(coords[1])])
                    fraccoords.append([-float(coords[0]), float(coords[1])])
                    fraccoords.append([float(coords[0]), -float(coords[1])])
                elif wp == "h" and coords is not None:
                    fraccoords.append([.5, float(coords[0])])
                    fraccoords.append([.5, -float(coords[0])])
                elif wp == "g" and coords is not None:
                    fraccoords.append([0, float(coords[0])])
                    fraccoords.append([0, -float(coords[0])])
                elif wp == "f" and coords is not None:
                    fraccoords.append([float(coords[0]), .5])
                    fraccoords.append([-float(coords[0]), .5])
                elif wp == "e" and coords is not None:
                    fraccoords.append([float(coords[0]), 0])
                    fraccoords.append([-float(coords[0]), 0])
                elif wp == "d":
                    fraccoords.append([.5, .5])
                elif wp == "c":
                    fraccoords.append([.5, 0])
                elif wp == "b":
                    fraccoords.append([0, .5])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == "pmg":
            for wp, coords in wyckoff_dict.items():
                if wp == "d" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[0]), -float(coords[1])])
                    fraccoords.append([-float(coords[0]) + .5, float(coords[1])])
                    fraccoords.append([float(coords[0]) + .5, -float(coords[1])])
                elif wp == "c" and coords is not None:
                    fraccoords.append([.25, float(coords[0])])
                    fraccoords.append([.75, -float(coords[0])])
                elif wp == "b":
                    fraccoords.append([0, .5])
                    fraccoords.append([.5, .5])
                elif wp == "a":
                    fraccoords.append([0, 0])
                    fraccoords.append([.5, 0])
        elif planegroup == "pgg":
            for wp, coords in wyckoff_dict.items():
                if wp == "c" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[0]), -float(coords[1])])
                    fraccoords.append([-float(coords[0]) + .5, float(coords[1]) + .5])
                    fraccoords.append([float(coords[0]) + .5, -float(coords[1]) + .5])
                elif wp == "b":
                    fraccoords.append([0, .5])
                    fraccoords.append([.5, 0])
                elif wp == "a":
                    fraccoords.append([0, 0])
                    fraccoords.append([.5, .5])
        elif planegroup == "cmm":
            for wp, coords in wyckoff_dict.items():
                if wp == "f" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[0]), -float(coords[1])])
                    fraccoords.append([-float(coords[0]), float(coords[1])])
                    fraccoords.append([float(coords[0]), -float(coords[1])])
                    fraccoords.append([float(coords[0]) + .5, float(coords[1]) + .5])
                    fraccoords.append([-float(coords[0]) + .5, -float(coords[1]) + .5])
                    fraccoords.append([-float(coords[0]) + .5, float(coords[1]) + .5])
                    fraccoords.append([float(coords[0]) + .5, -float(coords[1]) + .5])
                elif wp == "e" and coords is not None:
                    fraccoords.append([0, float(coords[0])])
                    fraccoords.append([0, -float(coords[0])])
                    fraccoords.append([.5, float(coords[0]) + .5])
                    fraccoords.append([.5, -float(coords[0]) + .5])
                elif wp == "d" and coords is not None:
                    fraccoords.append([float(coords[0]), 0])
                    fraccoords.append([-float(coords[0]), 0])
                    fraccoords.append([float(coords[0]) + .5, .5])
                    fraccoords.append([-float(coords[0]) + .5, .5])
                elif wp == "c":
                    fraccoords.append([.25, .25])
                    fraccoords.append([.75, .25])
                    fraccoords.append([.75, .75])
                    fraccoords.append([1.25, .75])
                elif wp == "b":
                    fraccoords.append([0, .5])
                    fraccoords.append([.5, 1])
                elif wp == "a":
                    fraccoords.append([0, 0])
                    fraccoords.append([.5, .5])
        elif planegroup == "p4":
            for wp, coords in wyckoff_dict.items():
                if wp == "d" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[0]), -float(coords[1])])
                    fraccoords.append([-float(coords[1]), float(coords[0])])
                    fraccoords.append([float(coords[1]), -float(coords[0])])
                elif wp == "c":
                    fraccoords.append([.5, 0])
                    fraccoords.append([0, .5])
                elif wp == "b":
                    fraccoords.append([.5, .5])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == "p4m":
            for wp, coords in wyckoff_dict.items():
                if wp == "g" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[0]), -float(coords[1])])
                    fraccoords.append([-float(coords[1]), float(coords[0])])
                    fraccoords.append([float(coords[1]), -float(coords[0])])
                    fraccoords.append([-float(coords[0]), float(coords[1])])
                    fraccoords.append([float(coords[0]), -float(coords[1])])
                    fraccoords.append([float(coords[1]), float(coords[0])])
                    fraccoords.append([-float(coords[1]), -float(coords[0])])
                elif wp == "f" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[0])])
                    fraccoords.append([-float(coords[0]), -float(coords[0])])
                    fraccoords.append([-float(coords[0]), float(coords[0])])
                    fraccoords.append([float(coords[0]), -float(coords[0])])
                elif wp == "e" and coords is not None:
                    fraccoords.append([float(coords[0]), .5])
                    fraccoords.append([-float(coords[0]), .5])
                    fraccoords.append([.5, float(coords[0])])
                    fraccoords.append([.5, -float(coords[0])])
                elif wp == "d" and coords is not None:
                    fraccoords.append([float(coords[0]), 0])
                    fraccoords.append([-float(coords[0]), 0])
                    fraccoords.append([0, float(coords[0])])
                    fraccoords.append([0, -float(coords[0])])
                elif wp == "c":
                    fraccoords.append([.5, 0])
                    fraccoords.append([0, .5])
                elif wp == "b":
                    fraccoords.append([.5, .5])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == "p4g":
            for wp, coords in wyckoff_dict.items():
                if wp == "d" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[0]), -float(coords[1])])
                    fraccoords.append([-float(coords[1]), float(coords[0])])
                    fraccoords.append([float(coords[1]), -float(coords[0])])
                    fraccoords.append([-float(coords[0]) + .5, float(coords[1]) + .5])
                    fraccoords.append([float(coords[0]) + .5, -float(coords[1]) + .5])
                    fraccoords.append([float(coords[1]) + .5, float(coords[0]) + .5])
                    fraccoords.append([-float(coords[1]) + .5, -float(coords[0]) + .5])
                elif wp == "c" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[0]) + .5])
                    fraccoords.append([-float(coords[0]), -float(coords[0]) + .5])
                    fraccoords.append([-float(coords[0]) + .5, float(coords[0])])
                    fraccoords.append([float(coords[0]) + .5, -float(coords[0])])
                elif wp == "b":
                    fraccoords.append([0, .5])
                    fraccoords.append([.5, 0])
                elif wp == "a":
                    fraccoords.append([0, 0])
                    fraccoords.append([.5, .5])
        elif planegroup == "p3":
            for wp, coords in wyckoff_dict.items():
                if wp == "d" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[1]), float(coords[0]) - float(coords[1])])
                    fraccoords.append([-float(coords[0]) + float(coords[1]), -float(coords[0])])
                elif wp == "c":
                    fraccoords.append([2/3, 1/3])
                elif wp == "b":
                    fraccoords.append([1/3, 2/3])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == "p3m":
            for wp, coords in wyckoff_dict.items():
                if wp == "e" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[1]), float(coords[0]) - float(coords[1])])
                    fraccoords.append([-float(coords[0]) + float(coords[1]), -float(coords[0])])
                    fraccoords.append([-float(coords[1]), -float(coords[0])])
                    fraccoords.append([-float(coords[0]) + float(coords[1]), float(coords[1])])
                    fraccoords.append([float(coords[0]), float(coords[0]) - float(coords[1])])
                elif wp == "d" and coords is not None:
                    fraccoords.append([float(coords[0]), -float(coords[0])])
                    fraccoords.append([float(coords[0]), 2 * float(coords[0])])
                    fraccoords.append([-2 * float(coords[0]), -float(coords[0])])
                elif wp == "c":
                    fraccoords.append([2/3, 1/3])
                elif wp == "b":
                    fraccoords.append([1/3, 2/3])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == "p3m1":
            for wp, coords in wyckoff_dict.items():
                if wp == "d" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[1]), float(coords[0]) - float(coords[1])])
                    fraccoords.append([-float(coords[0]) + float(coords[1]), -float(coords[0])])
                    fraccoords.append([float(coords[1]), float(coords[0])])
                    fraccoords.append([float(coords[0]) - float(coords[1]), -float(coords[1])])
                    fraccoords.append([-float(coords[0]), -float(coords[0]) + float(coords[1])])
                elif wp == "c" and coords is not None:
                    fraccoords.append([float(coords[0]), 0])
                    fraccoords.append([0, float(coords[0])])
                    fraccoords.append([-float(coords[0]), -float(coords[0])])
                elif wp == "b":
                    fraccoords.append([1/3, 2/3])
                    fraccoords.append([2/3, 1/3])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == "p6":
            for wp, coords in wyckoff_dict.items():
                if wp == "d" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[1]), float(coords[0]) - float(coords[1])])
                    fraccoords.append([-float(coords[0]) + float(coords[1]), -float(coords[0])])
                    fraccoords.append([-float(coords[0]), -float(coords[1])])
                    fraccoords.append([float(coords[1]), -float(coords[0]) + float(coords[1])])
                    fraccoords.append([float(coords[0]) - float(coords[1]), float(coords[0])])
                elif wp == "c":
                    fraccoords.append([.5, 0])
                    fraccoords.append([0, .5])
                    fraccoords.append([.5, .5])
                elif wp == "b":
                    fraccoords.append([1/3, 2/3])
                    fraccoords.append([2/3, 1/3])
                elif wp == "a":
                    fraccoords.append([0, 0])
        elif planegroup == "p6m":
            for wp, coords in wyckoff_dict.items():
                if wp == "f" and coords is not None:
                    fraccoords.append([float(coords[0]), float(coords[1])])
                    fraccoords.append([-float(coords[1]), float(coords[0]) - float(coords[1])])
                    fraccoords.append([-float(coords[0]) + float(coords[1]), -float(coords[0])])
                    fraccoords.append([-float(coords[0]), -float(coords[1])])
                    fraccoords.append([float(coords[1]), -float(coords[0]) + float(coords[1])])
                    fraccoords.append([float(coords[0]) - float(coords[1]), float(coords[0])])
                    fraccoords.append([-float(coords[1]), -float(coords[0])])
                    fraccoords.append([-float(coords[0]) + float(coords[1]), float(coords[1])])
                    fraccoords.append([float(coords[0]), float(coords[0]) - float(coords[1])])
                    fraccoords.append([float(coords[1]), float(coords[0])])
                    fraccoords.append([float(coords[0]) - float(coords[1]), -float(coords[1])])
                    fraccoords.append([-float(coords[0]), -float(coords[0]) + float(coords[1])])
                elif wp == "e" and coords is not None:
                    fraccoords.append([float(coords[0]), -float(coords[0])])
                    fraccoords.append([float(coords[0]), 2 * float(coords[0])])
                    fraccoords.append([-2 * float(coords[0]), -float(coords[0])])
                    fraccoords.append([-float(coords[0]), float(coords[0])])
                    fraccoords.append([-float(coords[0]), -2 * float(coords[0])])
                    fraccoords.append([2 * float(coords[0]), float(coords[0])])
                elif wp == "d" and coords is not None:
                    fraccoords.append([float(coords[0]), 0])
                    fraccoords.append([0, float(coords[0])])
                    fraccoords.append([-float(coords[0]), -float(coords[0])])
                    fraccoords.append([-float(coords[0]), 0])
                    fraccoords.append([0, -float(coords[0])])
                    fraccoords.append([float(coords[0]), float(coords[0])])
                elif wp == "c":
                    fraccoords.append([.5, 0])
                    fraccoords.append([0, .5])
                    fraccoords.append([.5, .5])
                elif wp == "b":
                    fraccoords.append([1/3, 2/3])
                    fraccoords.append([2/3, 1/3])
                elif wp == "a":
                    fraccoords.append([0, 0])
                elif wp == "honeycomb":
                    fraccoords.append([1/3,1/3])
                    fraccoords.append([2/3,2/3])
        return np.array(fraccoords)
    print("Planegroup not found")
    return np.array([])

def generallatticevectors(planegroup: str) -> np.ndarray:
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
        np.array([[.93969,.34202,0], [.98481,.17365,0], [0,0,1]]),  # p1
        np.array([[1.5,0,0], [0,1,0], [0,0,1]]),                     # p2,pm,pg,pgg,pmm,pmg
        np.array([[.75,.5,0], [.75,-.5,0], [0,0,1]]),               # cm,cmm
        np.array([[1,0,0], [0,1,0], [0,0,1]]),                       # p4,p4m,p4g
        np.array([[1,0,0], [.5,sqrt(3)/2,0], [0,0,1]])              # p3,p3m,p3m1,p6,p6m
    ]
    


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
    
    print("Planegroup not found")
    return latticevectorsset[0]


if __name__ == "__main__":
    #all the testing
    pts = builda2dcrystal("p6", ["a","c","d"])#fractional coords
    latticevectors = generallatticevectors("p6")#lattice vectors

    cartpts = np.zeros_like(pts)#turns form fractional into cartesian coords
    for i in range(pts.shape[0]):
        cartpts[i,:] = latticevectors[0,:2]*pts[i,0] + latticevectors[1,:2]*pts[i,1]
    print(cartpts)
    twodstructureplotter(latticevectors,cartpts,8,latticeshown=False,speciesshown=False,bondshown=True)#plots the structure