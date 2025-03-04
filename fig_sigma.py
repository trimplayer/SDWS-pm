import numpy as np
from matplotlib import pyplot as plt

#points1 = np.genfromtxt("sigma_001.csv")
#points2 = np.genfromtxt("sigma_005.csv")
#points3 = np.genfromtxt("stt_sigma_e001.csv",  delimiter=",")
#points4 = np.genfromtxt("stt_sigma_e005.csv",  delimiter=",")

points1 = np.genfromtxt("sigma_001_negative_nlgoff.csv")

plt.plot(points1[:,0]*1e9, points1[:,1]*1e9,label='ansys 001')
#plt.plot(points2[:,0]*1e9, points2[:,1]*1e9,label='ansys 005')
#plt.plot(points3[0,:]*0.1, points3[1,:]/1e9,label='001 1')
#plt.plot(points3[0,:]*0.1, points3[2,:]/1e9,label='001 15')
#plt.plot(points4[0,:]*0.1, points4[1,:]/1e9,label='005 1')
#plt.plot(points4[0,:]*0.1, points4[2,:]/1e9,label='005 15')

plt.legend(loc='best')

plt.show()