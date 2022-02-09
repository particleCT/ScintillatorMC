from ROOT import *
import sys
import numpy as np
import matplotlib.pyplot as plt
def GetHistArray(hist):
    n     = [hist.GetNbinsZ()+2, hist.GetNbinsY()+2,hist.GetNbinsX()+2]
    dx    = hist.GetXaxis().GetBinWidth(0)
    n2    = (n[0])*(n[1])*(n[2])

    d     = hist.GetArray()
    d.SetSize(n2)
    histA = np.array(d)
    histA = np.reshape(histA,n)[1:-1,1:-1,1:-1]
    return histA,dx

f = TFile(sys.argv[1])
print f.ls()
Edep    = f.Get("Edep_Tot")
LET     = f.Get("LET_Tot")
Entries = f.Get("Entries_Tot")

Edep,dx    = GetHistArray(Edep)
Entries,dx = GetHistArray(Entries)
LET, dx    = GetHistArray(LET)

plt.imshow(Edep[150])
plt.show()

Edep_line    = np.mean(np.mean(Edep,axis=0),axis=0)
Entries_line = np.mean(np.mean(Entries,axis=0),axis=0)
LET_line = np.mean(np.mean(LET,axis=0),axis=0)
plt.plot(Edep_line)#Entries[150,150])
plt.plot(Entries_line)#/Entries[150,150])
plt.plot(LET_line)#/Entries[150,150])
plt.show()
