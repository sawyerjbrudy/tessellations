from crystalgenerator import generallatticevectors, builda2dcrystal_batch
from crystalplotter import twodstructureplotter_batch
from math import pi, sqrt, cos, sin
import numpy as np
from matplotlib import pyplot as plt
import itertools

planegroup = "p6m"
lats = generallatticevectors(planegroup)
potential_wyckoffs = ['a','b','c','d','e','f']
wyckoffcombos = []

for r in range(1,len(potential_wyckoffs)+1):
    for combo in itertools.combinations(potential_wyckoffs,r):
        wyckoffcombos.append(combo)



for combo in wyckoffcombos:
    wyckoffdict = {}


    i = 0
    necessary_repeats = 1
    dimval = 0
    d_val = [0.0]
    e_val = [0.0]
    f_val = [0.0,0.0]

    if 'd' in combo:
        if i == 0:
            necessary_repeats = necessary_repeats*4
            dimval+=1

    if 'e' in combo:
        if i == 0:
            necessary_repeats = necessary_repeats*4
            dimval += 1
        
    if 'f' in combo:
        if i == 0:
            necessary_repeats = necessary_repeats*16
            dimval += 2

    if dimval > 0:
        valarray = np.empty([4]*dimval, dtype=object)
        for idx, val in np.ndenumerate(valarray):
            valarray[idx] = list((np.array(idx)+1)*2/10)
        flatvalarray = valarray.flatten()

    while i < necessary_repeats:
        if 'a' in combo:
            wyckoffdict["a"] = None
        if 'b' in combo:
            wyckoffdict["b"] = None
        if 'c' in combo:
            wyckoffdict["c"] = None
        
        if dimval > 0:
            coord_list = flatvalarray[i]  # This is a list with all needed coordinates
            
            # Assign coordinates to present Wyckoffs in order
            coord_index = 0
            if 'd' in combo:
                wyckoffdict["d"] = [coord_list[coord_index]]
                coord_index += 1
            if 'e' in combo:
                wyckoffdict["e"] = [coord_list[coord_index]]
                coord_index += 1
            if 'f' in combo:
                wyckoffdict["f"] = coord_list[coord_index:coord_index+2]  # f needs 2 coordinates
                coord_index += 2
        


        fracs = builda2dcrystal_batch(planegroup,wyckoff_dict=wyckoffdict,)
        twodstructureplotter_batch(lats, fracs,4,speciesshown=True,bondshown=True,plane=planegroup, wyckoffdict=wyckoffdict, folder="p6mpossibilities")

        i += 1