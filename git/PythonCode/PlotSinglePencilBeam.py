from root_numpy import tree2array
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import ROOT
from ROOT import *
import sys
import matplotlib.pyplot as plt

rootfile   = TFile(sys.argv[1])
Hist3D     = rootfile.Get("Edep_Tot")
phasespace = tree2array(rootfile.Get("PencilBeam"), branches=["idPB","binglobal","Estop"])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ids = np.where(phasespace["idPB"]==int(sys.argv[2]))
phasespace = phasespace[ids]
bx,by,bz = (ROOT.Long(),ROOT.Long(),ROOT.Long())
x =[]
y =[]
z =[]
E =[]
for idPB, binglobal, Estop in phasespace:
    Hist3D.GetBinXYZ(int(binglobal),bx,by,bz)
    x.append(Hist3D.GetXaxis().GetBinCenter(bx))
    y.append(Hist3D.GetYaxis().GetBinCenter(by))
    z.append(Hist3D.GetZaxis().GetBinCenter(bz))
    E.append(Estop)
    
sc = ax.scatter(x,y,z,c=E)
ax.set_xlim([-150,150])
ax.set_ylim([-150,150])
ax.set_zlim([-150,150])
plt.colorbar(sc)
plt.show()
