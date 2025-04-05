import numpy as np
from matplotlib import pyplot as plt


points1= np.genfromtxt("surface_y2.csv", delimiter=",")
points2= np.genfromtxt("surface_y15.csv", delimiter=",")
points3= np.genfromtxt("surface_y07.csv", delimiter=",")
points4= np.genfromtxt("surface_y05.csv", delimiter=",")


plt.plot(points1[0,:],points1[1,:],label='y = 2.0',color='black')
plt.plot(points2[0,:],points2[1,:],label='y = 1.5',color='black')
plt.plot(points3[0,:],points3[1,:],label='y = 0.7',color='black')
plt.plot(points4[0,:],points4[1,:],label='y = 0.5',color='black')



plt.grid(True)
plt.xlabel(r"$\frac{{{x_1}}}{a}$",fontsize=16)
plt.ylabel(r"$\frac{{\varepsilon f}}{a}$",fontsize=16,rotation=0)

# Set custom y-ticks and labels
custom_yticks = [-1, -0.5, 0, 0.5, 1]  # Your desired tick positions
custom_yticklabels = [r"$- \varepsilon$", r"$-0.5 \varepsilon$", "0", r"$0.5 \varepsilon$", r"$ \varepsilon$"]  # Corresponding labels

plt.yticks(custom_yticks, custom_yticklabels)

# Add text annotations
plt.text(0.5, 0.97, r'$1$', fontsize=10, color='black')
plt.text(0.48, 0.8, r'$2$', fontsize=10, color='black')
plt.text(0.49, 0.65, r'$3$', fontsize=10, color='black')
plt.text(0.49, 0.4, r'$4$', fontsize=10, color='black')

plt.savefig("surface_plot.png", dpi=300, bbox_inches='tight')

plt.show()