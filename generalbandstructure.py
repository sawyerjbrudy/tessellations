import numpy as np
import matplotlib.pyplot as plt
import pythtb as tb
from crystalgenerator import builda2dcrystal,generallatticevectors
from crystalplotter import twodstructureplotter, supercellmaker
from nearestneighbor import nearestneighborfinder
from seekpath import get_path
from hoppingparameterfinder import hopindexer



def bandstructure(plane,wyckoffs,hop_param):

  lat = generallatticevectors(plane)
  orb = builda2dcrystal(plane,wyckoffs)
  lat2d = lat[:2,:2]


  onsite_param = 1

  inds = []
  for i in range(orb.shape[0]):
      inds.append(i)
  tbm = tb.tb_model(2,2,lat2d,orb)
  tbm.set_onsite([onsite_param]*orb.shape[0])

  #hopping parameters
  hoppingindex = hopindexer(orb,lat)
  for i in hoppingindex:
      tbm.set_hop(hop_param,i[0],i[1],[i[2][0],i[2][1]],allow_conjugate_pair=True)

  #using seekpath(must cite) to get the kpath
  pos = np.column_stack((orb,np.zeros(orb.shape[0])))
  species = [6]*orb.shape[0]
  structure = (lat,pos,species)

  pathunformatted = get_path(structure)['point_coords']
  path = []
  labels = []
  for key,value in pathunformatted.items():
      if value[0:2] not in path:
          path.append(value[0:2])
          labels.append(key)

  #complete the path
  path.append(path[0])
  labels.append(labels[0])

  print(path)
  print(labels)


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

  ax.set_title("Bandstructure of "+plane+" "+str(wyckoffs)+" with hop param = "+str(hop_param))
  ax.set_xlabel("k")
  ax.set_ylabel("E")

  plt.show()





  # tbm.display()

  # # visualize the model
  # c1 = tbm.cut_piece(8,0)
  # c2 = c1.cut_piece(8,1)

  # fig = c2.visualize(0,1)
  # plt.show()

if __name__ == "__main__":
    bandstructure("p6m",["c"],-.5)