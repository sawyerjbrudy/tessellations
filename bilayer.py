from crystalgenerator import generallatticevectors, builda2dcrystal
from crystalplotter import twodstructureplotter, supercellmaker
from math import pi, sqrt, cos, sin
import numpy as np
from matplotlib import pyplot as plt

import itertools

def moireprimitivecell(latticevectors1, latticevectors2, misorientationangle):
    R = 10
    t = .05
    twodlatA = np.array([latticevectors1[0,:2],latticevectors1[1,:2]])
    twodlatB = np.array([latticevectors2[0,:2],latticevectors2[1,:2]])

    Aexact = np.linalg.solve(twodlatA, twodlatB)
    

    I = []
    for vals in itertools.product(range(-R,R+1), repeat=8):
        M1 = np.array([[vals[0],vals[1]],[vals[2],vals[3]]])
        M2 = np.array([[vals[4],vals[5]],[vals[6],vals[7]]])
        try:
            M1i = np.linalg.inv(M1)
            Apot = M1i @ M2
            if np.max(np.abs(Apot - Aexact)) < t:
                I.append(M2)
                break
        except np.linalg.LinAlgError:
            continue

    print("found solutions within tolerance")

    sopts = np.zeros(len(I)) 
    i = 0
    for n in I:
        dets = np.linalg.det(n)
        Alats = np.linalg.norm(np.cross(twodlatB[0],twodlatB[1]))

        sopts[i] = dets*Alats
        i += 1
    
    sminind = np.where(sopts == np.min(sopts))[0]

    minimizedI = I[sminind]
    primlats = minimizedI @ twodlatB


    return primlats

def moireanycell(latticevectors1, latticevectors2, misorientationangle):
    R = 10
    t = .05
    twodlatA = np.array([latticevectors1[0,:2],latticevectors1[1,:2]])
    twodlatB = np.array([latticevectors2[0,:2],latticevectors2[1,:2]])

    Aexact = np.linalg.solve(twodlatA, twodlatB)
    

    I = []
    for vals in itertools.product(range(-R,R+1), repeat=8):
        M1 = np.array([[vals[0],vals[1]],[vals[2],vals[3]]])
        M2 = np.array([[vals[4],vals[5]],[vals[6],vals[7]]])
        try:
            M1i = np.linalg.inv(M1)
            Apot = M1i @ M2
            if np.max(np.abs(Apot - Aexact)) < t:
                I.append(M2)
                break
        except np.linalg.LinAlgError:
            continue

    print("found solutions within tolerance")

    sopts = np.zeros([len(I),1]) 
    i = 0
    for n in I:
        dets = np.linalg.det(n)
        Alats = np.linalg.norm(np.cross(twodlatB[0],twodlatB[1]))

        sopts[i] = dets*Alats
        i += 1
    
    sminind = np.where(sopts == np.min(sopts))[0]


    minimizedI = I[sminind]
    primlats = minimizedI @ twodlatB


    return primlats

moireprimitivecell(generallatticevectors("p6m"),generallatticevectors("p6m"),0)



