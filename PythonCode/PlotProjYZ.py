from ROOT import *
import sys
import numpy as np
import matplotlib.pyplot as plt
def GetHistArray(hist):
    n     = [hist.GetNbinsY()+2,hist.GetNbinsX()+2]
    dx    = hist.GetXaxis().GetBinWidth(0)
    n2    = (n[0])*(n[1])

    d     = hist.GetArray()
    d.SetSize(n2)
    histA = np.array(d)
    histA = np.reshape(histA,n)[1:-1,1:-1]
    return histA,dx

f = TFile(sys.argv[1])
print f.ls()
A = np.zeros((10000, 300, 300),dtype=np.float32)
for i in range(0,100):
    if(i%100==0): print i
    Edep       = f.Get("YZProj_Q/YZProj_Q_"+str(i))
    Edep,dx    = GetHistArray(Edep)
    A[i] = Edep
    #np.savetxt("SlantedEdge180Proj/YZProj_"+str(i)+".txt",Edep)
np.savez("SlantedEdge180.npz",A)
