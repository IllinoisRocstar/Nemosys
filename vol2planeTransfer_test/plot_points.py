from scipy import loadtxt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

x1, y1, rho = loadtxt("has_sphere_has_coord", unpack=True, skiprows=2, usecols=(1,2,4))
fig = plt.figure()
ax = fig.gca(projection='3d')
x1_zeroes = []
x1_non_zeroes = []
y1_zeroes = []
y1_non_zeroes = []
#for i in range(0,len(x1[0:2822])):
#	if (rho[i] == 0) :
#		x1_zeroes.append(x1[i])
#		y1_zeroes.append(y1[i])
for i in range(0,len(x1)):
		if (rho[i] == 0):
			x1_zeroes.append(x1[i])
			y1_zeroes.append(y1[i])
		else:
			x1_non_zeroes.append(x1[i])
			y1_non_zeroes.append(y1[i])


#x1_all_zeroes = x1[2823:-1]
#y1_all_zeroes = y1[2823:-1]
#print len(x1_zeroes)
#print len(x1_non_zeroes)
#print len(x1_all_zeroes)
ax.plot(x1_zeroes,y1_zeroes,.8, 'b.')
ax.plot(x1_non_zeroes, y1_non_zeroes,.8, 'r.')
#ax.plot(x1_all_zeroes, y1_all_zeroes,.8, 'g.')
plt.show()
#ax.plot(x1, -1*fullUx*gradphi[i][:], line_fmts2[int(floor(i/10))], linewidth = 1, label=label2) 	
#ax.set_xlabel('X ($m$)',fontsize=15)
#ax.set_ylabel('$\\phi$',fontsize=15)
#ax.grid(True)

#plt.legend(loc = 'lower left', prop={'size':13}, ncol=2, fancybox=False, shadow=False, frameon=True,markerscale=1.5)
#plt.show()
