from scipy import loadtxt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

x1, y1, rho = loadtxt("plane_interp.out", unpack=True, skiprows=2, usecols=(1,2,4))
fig = plt.figure()
ax = fig.gca(projection='3d')
x1_zeroes = []
x1_others = []
x1_non_zeroes = []
y1_zeroes = []
y1_others = []
y1_non_zeroes = []

for i in range(0,len(x1)):
		if (rho[i] == 0 and i < len(x1) - 7):
			x1_zeroes.append(x1[i])
			y1_zeroes.append(y1[i])
		elif (rho[i] == 0 and i >= len(x1) - 7):
			x1_others.append(x1[i])
			y1_others.append(y1[i])
		else:
			x1_non_zeroes.append(x1[i])
			y1_non_zeroes.append(y1[i])


ax.plot(x1_zeroes,y1_zeroes,.8, 'b.')
ax.plot(x1_non_zeroes, y1_non_zeroes,.8, 'r.')
ax.plot(x1_others, y1_others, .8, 'g.')
plt.show()
