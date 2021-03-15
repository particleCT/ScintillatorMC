from ROOT import *
import numpy as np
import sys
from root_numpy import hist2array
f=TFile(sys.argv[1])
print(f.ls())
hist = f.Get("YZProj_Tot_Q")
hist = hist2array(hist)
np.savetxt("YZProj_Tot_Q.txt",hist)
