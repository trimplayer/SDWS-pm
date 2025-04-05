import numpy as np
from matplotlib import pyplot as plt


points1= np.genfromtxt("stt_fou_y2.csv", delimiter=",")
points2= np.genfromtxt("stt_fou_y15.csv", delimiter=",")
points3= np.genfromtxt("stt_fou_y07.csv", delimiter=",")
points4= np.genfromtxt("stt_fou_y05.csv", delimiter=",")



plt.plot(points1[0,:],points1[1,:]/1e9,label=r"$1- y=2.0$", color='black')
plt.plot(points2[0,:],points2[1,:]/1e9,label=r"$2- y=1.5$", color='black')
plt.plot(points3[0,:],points3[1,:]/1e9,label=r"$3- y=0.7$", color='black')
plt.plot(points4[0,:],points4[1,:]/1e9,label=r"$4- y=0.5$", color='black')

# Add text annotations
plt.text(0.1, 0.041, r'$1$', fontsize=10, color='black')
plt.text(0.098, -0.146, r'$2$', fontsize=10, color='black')
plt.text(0.099, -0.88, r'$3$', fontsize=10, color='black')
plt.text(0.1, -1.95, r'$4$', fontsize=10, color='black')


leg = plt.legend(loc='best', handlelength=0, handletextpad=0, fancybox=True,fontsize=14)
for item in leg.legendHandles:
    item.set_visible(False)
plt.xlabel(r"$\varepsilon$",fontsize = 16)
plt.ylabel(r"${\sigma _{tt}}$, GPa",fontsize = 16)
plt.grid()
plt.savefig("fou_fig.png", dpi=300, bbox_inches='tight')
plt.show()
