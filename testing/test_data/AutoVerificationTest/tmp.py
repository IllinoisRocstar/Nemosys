from __future__ import print_function
from pyNemosys import *
import numpy as np
import matplotlib.pyplot as plt

def errorfit(m,x,b):
  return m*x + b

fine = meshBase.Create('fine.vtu')
coarse = meshBase.Create('coarse.vtu')
finer = meshBase.Create('finer.vtu')
#finer = meshBase.Create('f1.vtu')
#fine = meshBase.Create('f2.vtu')
#coarse = meshBase.Create('f3.vtu')


#x = []
#x0 = sum(finer.getCellLengths())/finer.getNumberOfCells()
#x1 = sum(fine.getCellLengths())/fine.getNumberOfCells()
#x2 = sum(coarse.getCellLengths())/coarse.getNumberOfCells()
#x.append(np.log(x0))
#x.append(np.log(x1))
#x.append(np.log(x2))
#print(x)

arr = intV(1)
arr[0] = 0

oac = OrderOfAccuracy(finer,fine,coarse,arr)
oac.computeRichardsonExtrapolation()
finer.write("richardson.vtu")
#oac.computeMeshWithResolution(.005, "richierich.vtu")

resolution = oac.computeResolution(.005)
gci_21 = doubleVV()
gci_32 = doubleVV()
gci_21 = oac.computeGCI_21()
gci_32 = oac.computeGCI_32()

orderOfAcc = doubleVV()
orderOfAcc = oac.getOrderOfAccuracy()

asymp = doubleVV()
asymp = oac.checkAsymptoticRange()

def p(cap, r, c, x):
  print(cap)
  for j in range(0,r):
    for k in range(0,c):
      print(x[j][k])

p('gci_21',1,3,gci_21)
p('gci_32',1,3,gci_32)
p('asymp ratios',1,3,asymp)
p('oac',1,3,orderOfAcc)
p('resolution',1,3,resolution)
