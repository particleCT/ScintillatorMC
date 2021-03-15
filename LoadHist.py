from ROOT import *
from root_numpy import hist2array
import sys
import numpy as np
import matplotlib.pyplot as plt

f=TFile(sys.argv[1])
hist = f.Get("YXProj/YXProj_0")
print(hist.GetXaxis().GetBinCenter(76))
print(hist.GetXaxis().GetBinLowEdge(76))
print(hist.GetXaxis().GetBinUpEdge(76))
hist = hist2array(hist)
print(hist.shape)
np.save("test.npy",hist)
