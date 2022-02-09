from ROOT import *
import numpy as np
import sys
def GetHistArray(hist):
    n     = [hist.GetNbinsZ()+2,hist.GetNbinsY()+2,hist.GetNbinsX()+2]
    dx    = hist.GetXaxis().GetBinWidth(0)
    n2    = (n[0])*(n[1])*(n[2])
    
    d     = hist.GetArray()
    d.SetSize(n2)
    histA = np.array(d)
    histA = np.reshape(histA,n)[1:-1,1:-1,1:-1]
    return histA
                    
f = TFile(sys.argv[1])
hist = f.Get("L_Tot")
hist = GetHistArray(hist)
print hist.shape
np.savetxt("Reference200MeV.txt",hist.flatten())
