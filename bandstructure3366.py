import numpy as np
import matplotlib.pyplot as plt
import pythtb as tb
from crystalgenerator import builda2dcrystal,generallatticevectors
from crystalplotter import twodstructureplotter
from seekpath import get_path
from pymatgen.core import Structure, Lattice

"""
Band Structure Calculation for 3366 Tessellation Pattern

This script calculates and visualizes the electronic band structure for a 2D crystal
with p6m symmetry and Wyckoff position "c", implementing a tight-binding model
with nearest-neighbor hopping interactions.

The code constructs a tight-binding model for a hexagonal lattice structure
(3366 tessellation) and computes the energy eigenvalues along a predefined
k-path in the Brillouin zone, then plots the resulting band structure.

Dependencies
-----------
numpy : Numerical computing library
matplotlib.pyplot : Plotting library
pythtb : Tight-binding model library
tssltnwyckoff : Custom module for 2D crystal construction
WyckoffPlaneGroupID : Custom module for structure plotting
seekpath : K-path generation (imported but not used)
pymatgen.core : Structure and lattice handling (imported but not used)

Parameters
----------
hop_param : float, default=0.5
    Hopping parameter for nearest-neighbor interactions in the tight-binding model.
    Controls the strength of electron hopping between adjacent sites.

onsite_param : float, default=0
    On-site energy parameter for all atomic orbitals. Set to zero for simplicity.

Crystal Structure
----------------
- Space group: p6m (plane group)
- Wyckoff position: "c"
- Lattice: 2D hexagonal lattice
- Number of orbitals: Determined by builda2dcrystal function

Tight-Binding Model
-------------------
- Dimension: 2D
- Number of orbitals: Variable (determined by crystal structure)
- Hopping terms: Nearest-neighbor interactions with specified hop_param
- Hopping connections:
  * Site 0 ↔ Site 1 (same unit cell)
  * Site 1 ↔ Site 2 (same unit cell) 
  * Site 0 ↔ Site 2 (same unit cell)
  * Site 2 ↔ Site 0 (adjacent unit cell [0,1])
  * Site 2 ↔ Site 1 (adjacent unit cell [1,0])
  * Site 0 ↔ Site 1 (adjacent unit cell [1,-1])

K-Path
------
Manual k-path through high-symmetry points:
- Γ (0, 0): Gamma point
- M (1/3, 1/3): M point  
- K (1/2, 0): K point
- Γ (0, 0): Return to Gamma point
- Number of k-points: 301

Output
------
- Band structure plot showing energy eigenvalues vs. k-vector
- X-axis: k-vector distance with labeled high-symmetry points
- Y-axis: Energy eigenvalues
- Vertical lines: High-symmetry point markers
- Title: Includes hopping parameter value

Notes
-----
- The code includes commented sections for alternative k-path generation
  using seekpath library (requires citation)
- Visualization of the tight-binding model structure is available
  in commented code sections
- The "3366" in the title refers to the tessellation pattern type

References
----------
- seekpath library should be cited if used for k-path generation
- pythtb library for tight-binding calculations
- Custom modules for crystal structure generation
"""

orb = builda2dcrystal("p6m",["c"])
lat = generallatticevectors("p6m")
lat2d = lat[:2,:2]
print(lat)

hop_param = -.1
onsite_param = 5

inds = []
for i in range(orb.shape[0]):
    inds.append(i)
tbm = tb.tb_model(2,2,lat2d,orb)
tbm.set_onsite([onsite_param]*orb.shape[0])

#manual hopping parameters
tbm.set_hop(hop_param,0,1,[0,0])
tbm.set_hop(hop_param,1,2,[0,0])
tbm.set_hop(hop_param,0,2,[0,0])
tbm.set_hop(hop_param,2,0,[0,1])
tbm.set_hop(hop_param,2,1,[1,0])
tbm.set_hop(hop_param,0,1,[1,-1])

# #using seekpath(must cite) to get the kpath
# species = ['A']*orb.shape[0]
# structure = Structure(Lattice(lat), species, orb)
# kpath = get_path(structure)
# print(kpath)

#manual kpath
path = [[0,0],[.33333,.33333],[.5,0],[0,0]]
labels = (r"$\Gamma$",r"$M$",r"$K$",r"$\Gamma$")
(k_vec,k_dist,k_node) = tbm.k_path(path,301)

#plotting the bandstructure
evals = tbm.solve_all(k_vec)

# First make a figure object
fig, ax = plt.subplots()

# specify horizontal axis details
ax.set_xlim(k_node[0],k_node[-1])
ax.set_xticks(k_node)
ax.set_xticklabels(labels)
for n in range(len(k_node)):
  ax.axvline(x=k_node[n], linewidth=0.5, color='k')

# plot bands
for n in range(len(evals)):
  ax.plot(k_dist,evals[n])

ax.set_title("Bandstructure of 3366 with hop param = "+str(hop_param))
ax.set_xlabel("k")
ax.set_ylabel("E")

plt.show()





# tbm.display()

# #visualize the model
# c1 = tbm.cut_piece(8,0)
# c2 = c1.cut_piece(8,1)

# fig = c2.visualize(0,1)
# plt.show()

