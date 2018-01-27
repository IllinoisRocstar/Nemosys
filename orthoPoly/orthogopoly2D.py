from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import time

class OrthoPoly:

  def __init__(self, order, x, fun):
    self.order = order
    self.a = np.zeros(order+1)
    self.b = np.zeros(order)
    self.x = x
    self.y = np.zeros((order+1,len(x))) 
    self.f = fun(x)
    self.phiTphiInv = self.ComputePhiTPhiInv()
    
  def EvaluateOrthogonal(self, power, xk):
    p0=1.
    if power == 0:
      return p0
    p1 = xk - self.a[0]
    if power == 1:
      return p1
    for m in range(1,power):
      p2 = (xk - self.a[m])*p1 - self.b[m-1]*p0
      p0 = p1
      p1 = p2
    return p2 
 
  def ComputeAB(self): 
    self.a[0] = 0.
    n = len(self.x)
    for k in range(0,n):
      self.a[0] += self.x[k]
    self.a[0] /= n
    for m in range(1, self.order):
      sum0 = sum1 = sum2 = 0.
      for k in range(0,n):
        tmp0 = self.EvaluateOrthogonal(m-1,self.x[k])
        tmp1 = self.EvaluateOrthogonal(m,self.x[k])
        sum0 += tmp0*tmp0
        sum1 += tmp1*tmp1
        sum2 += self.x[k]*tmp1*tmp1
      self.a[m] = sum2/sum1
      self.b[m-1] = sum1/sum0
  
  def EvaluateOrthogonals(self):
    self.ComputeAB()
    for i in range(self.order+1):
      for j in range(len(self.x)):
        self.y[i][j] = self.EvaluateOrthogonal(i,x[j])
  
  def PlotPoly(self):
    #self.EvaluateOrthogonals()
    ax = plt.axes()   
    for i in range(self.order+1):
      ax.plot(self.x,self.y[i])
    plt.show()

  def ComputePhiTPhiInv(self):
    self.EvaluateOrthogonals()
    phi = np.transpose(self.y)
    phiTphi = np.matmul(self.y,phi)
    phiTphiDiag = np.diag(phiTphi)
    phiTphiDiagInv = [0 for i in range(len(phiTphiDiag))]
    for j in range(len(phiTphiDiag)):
      phiTphiDiagInv[j] = 1./phiTphiDiag[j]
    return np.diagflat(phiTphiDiagInv)
 
  def ComputeFit(self):
    a = np.matmul(np.matmul(self.phiTphiInv,self.y) ,f)
    self.fapprox = np.matmul(phi,a)
  
  def PlotFit(self):
    self.ComputeFit()
    ax = plt.axes()
    ax.plot(self.x, self.f, 'b')
    ax.plot(self.x, self.fapprox,'m.')
    plt.show()

def Compute3DPolyFit(power, x, y, z, FDummy, G):
    X, Y, Z = np.meshgrid(x,y,z)
    F = G(X,Y,Z)
    np.savetxt('F.txt',np.ndarray.flatten(F))
    start1 = time.time()
    if power == 1:
    	toremove = [3,5,6,7]
    elif power == 2:
    	toremove = [5,7,8,11,13,14,15,16,17,19,20,21,22,23,24,25,26]
    elif power == 3:
    	toremove = [7,10,11,13,14,15,19,22,23,25,26,27,28,29,30,31,34,35,37,38,39,40,41,42,43,44,45,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]
    
    opx = OrthoPoly(power,x,FDummy)
    opy = OrthoPoly(power,y,FDummy)
    opz = OrthoPoly(power,z,FDummy)

    kronProd = np.kron(np.matmul(opx.phiTphiInv,opx.y),np.matmul(opy.phiTphiInv,opy.y))
    kronProd = np.kron(kronProd, np.matmul(opz.phiTphiInv,opz.y))
    
    #print opx.y.shape[0], opx.y.shape[1],kronProd
    
    #print len(toremove)
    start = time.time()
    kronProd = np.delete(kronProd,toremove,0)
    print "kronProd shape: ", kronProd.shape[0], kronProd.shape[1]
    
    a = np.matmul(kronProd,np.ndarray.flatten(F))
    end = time.time()
    print "Time elapsed: ", (end-start)*1000.0
    print a
    #NewPhi = np.kron(np.transpose(opx.y),np.transpose(opy.y))
    #NewPhi = np.kron(NewPhi, np.transpose(opz.y))
    #NewPhi = np.delete(NewPhi,toremove,1)
    #FApprox = np.matmul(NewPhi,a)
    #Error = np.absolute(FApprox - np.ndarray.flatten(F))
    #
    result = []
    #result.append(FApprox)
    #result.append(Error)
    phix = np.zeros(power+1)
    phiy = np.zeros(power+1)
    phiz = np.zeros(power+1)
    for i in range(power+1):
      phix[i] = opx.EvaluateOrthogonal(i,.05)
      phiy[i] = opy.EvaluateOrthogonal(i,.5)  
      phiz[i] = opz.EvaluateOrthogonal(i,.5)  
    
    Phi = np.kron(phix, phiy)
    Phi = np.kron(Phi, phiz)
    Phi = np.delete(Phi,toremove,0) 
    print Phi
    print np.matmul(Phi,a)
    end1 = time.time()
    print "total time: ", (end1-start1)*1000.0
    return result 


def Compute2DPolyFit(power, x, y, FDummy, G):
  opx = OrthoPoly(power,x,FDummy)
  opy = OrthoPoly(power,y,FDummy)
  kronProd = np.kron(np.matmul(opx.phiTphiInv,opx.y),np.matmul(opy.phiTphiInv,opy.y))
  print kronProd.shape[0], kronProd.shape[1]
  #numRow = kronProd.shape[0]
  #numRemove = power*(power+1)/2
  #toremove = [5,7,8]
  #for k in range(numRemove):
  #  toremove.append(numRow - 1 - k)
  #print toremove
  #kronProd = np.delete(kronProd,toremove,0)
  print kronProd.shape[0], kronProd.shape[1]
  X, Y = np.meshgrid(x,y)
  Z = G(X,Y)
  a = np.matmul(kronProd,np.ndarray.flatten(Z))
  print a
  NewPhi = np.kron(np.transpose(opx.y),np.transpose(opy.y))
  #NewPhi = np.delete(NewPhi,toremove,1)
  print NewPhi.shape[0], NewPhi.shape[1]
  ZApprox = np.matmul(NewPhi,a)
  Error = np.absolute(ZApprox - np.ndarray.flatten(Z))
  result = []
  result.append(ZApprox)
  result.append(Error)
  return result 
  
def F(x):
 return np.divide((np.add(x,1./10))*np.sin(np.add(5*x,-1)),np.add(1,np.square(x*(np.sin(x-1/2.)))))

x = np.arange(-1,1,0.015)
y = np.arange(-1,1,0.015)
z = np.arange(-1,1,0.015)
#X, Y, Z = np.meshgrid(x,y,z)

#def G(X,Y,Z):
#  return 467.4*X + 128*Y + 251.1*Z + 27.*X*Y + 1823*X**2 + .723*Y**2 + .1938*X*Z + 670*Y*Z + 3000*Z**2 + .1782*X*Y*Z + 21839*X*Z**2 - 18234*Z*Y**2

#result = Compute3DPolyFit(3,x,y,z,F,G)
#print "max error ", max(result[1])
#print "min error ", min(result[1])


#def G(x,y):
#  return 467.4*X + .2*X**2 + 1.234*Y + .5*X*Y + .3*Y**2
#  #
def G(x,y):
	return np.multiply(F(x),F(y))
#
X, Y = np.meshgrid(x, y)
Z = G(X,Y) #F(X)*F(Y)
result = Compute2DPolyFit(15,x,y,F,G)
print "max error ", max(result[1])
print "min error ", min(result[1])
#
#
#
fig = plt.figure()
ax = fig.gca(projection='3d')
#
# Plot the surface.
#surf = ax.plot_surface(X, Y, np.reshape(result[0],(len(x),len(x))), cmap=cm.coolwarm,
#                       linewidth=0)#, antialiased=False)
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
#
#surf = ax.plot_surface(X, Y, np.reshape(result[1],(len(x),len(x))), cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)
#
#
# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#
## Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
#
plt.show()
