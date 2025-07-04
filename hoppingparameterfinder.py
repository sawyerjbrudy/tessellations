from crystalplotter import supercellmaker
from nearestneighbor import nearestneighborfinder
from crystalgenerator import builda2dcrystal,generallatticevectors
import numpy as np

def hopindexer(fraccoords,latticevectors):
   supercell = supercellmaker(latticevectors,fraccoords)
   neighborfinder = supercell.copy().reshape(-1,2)
   nearestneighbors = nearestneighborfinder(neighborfinder)
   hoppingindex = []
   for i in range(supercell.shape[2]):
      for j in nearestneighbors:

        if np.allclose(supercell[3,3,i], j[0:2]):
            # Find where j[2:4] appears in the supercell
            searchspec = j[2:4]
            found = np.where(np.all(np.isclose(supercell, searchspec, atol=1e-5), axis=-1))

            # Only proceed if matches are found
            if len(found[0]) > 0 and len(found[1]) > 0 and len(found[2]) > 0:
                cell = [found[0][0], found[1][0]]
                secondindex = found[2][0]
                truecell = np.array(cell)
                perfectcell = truecell - np.array([3,3])
                hoppingindex.append([i, secondindex, list(perfectcell)])
   return hoppingindex

if __name__ == "__main__":
   orb = builda2dcrystal("p6m",["c"])
   lat = generallatticevectors("p6m")
   print(hopindexer(orb,lat))

