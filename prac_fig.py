import numpy as np
from matplotlib import pyplot as plt

points1= np.genfromtxt("stt_a100_Ms60_1.csv", delimiter=",")
points2= np.genfromtxt("stt_a100_Ms60_T0.csv", delimiter=",")
points3= np.genfromtxt("stt_a100_Ms60_sig0.csv", delimiter=",")

#a  = 100

fig, ax = plt.subplots(nrows=1, ncols=2,layout="constrained")
ax[0].tick_params(axis='both', which='major', labelsize=18)
ax[1].tick_params(axis='both', which='major', labelsize=18)

ax[0].plot(points1[0, 4:10], points1[-1, 4:10] /1e9,'r',label="$T=0.1$ ГПа; ${\sigma _{0}^{s}}=1$ Н/м")
ax[0].plot(points2[0, 4:10], points2[-1, 4:10] /1e9,'g',label="$T=0$ ГПа; ${\sigma _{0}^{s}}=1$ Н/м" )
ax[0].plot(points3[0, 4:10], points3[-1, 4:10] /1e9, 'b',label="$T=0.1$ ГПа; ${\sigma _{0}^{s}}=0$ Н/м")

ax[0].plot(points1[0, 4:10], points1[1, 4:10] /1e9,'r--')
ax[0].plot(points2[0, 4:10], points2[1, 4:10] /1e9,'g--')
ax[0].plot(points3[0, 4:10], points3[1, 4:10] /1e9,'b--')

#ax[0].annotate(r"$1$", xy=(12.5,0.3), xytext=(14.5,0.25),arrowprops=dict(arrowstyle="-",facecolor='black'),fontsize=20)
#ax[0].annotate(r"$1$", xy=(12.5,0.197), xytext=(14.5,0.25),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)

#ax[0].annotate(r"$2$", xy=(11.5,0.2437), xytext=(12,0.2193),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)
#ax[0].annotate(r"$2$", xy=(10.5,0.1799), xytext=(12,0.2193),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)

#ax[0].annotate(r"$3$", xy=(14,0.056), xytext=(15,0.03),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)
#ax[0].annotate(r"$3$", xy=(14,0.0028), xytext=(15,0.03),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)

#leg = ax[0].legend(['1 - $T=0.1$ GPa; ${\sigma _{0}^{s}}=1$ N/m','2 - $T=0$ GPa; ${\sigma _{0}^{s}}=1$ N/m', '3 - $T=0.1$ GPa; ${\sigma _{0}^{s}}=0$ N/m'],loc='best',handlelength=0, handletextpad=0, fancybox=True)
#for item in leg.legendHandles:
#    item.set_visible(False)
#plt.title("snn a=100")


#xi = list(range(len(points1[0,4:15]*a)))
ax[0].set_xticks(points1[0,4:10:2], [ '%.2f' % elem for elem in points1[0,4:10:2]])

ax[0].grid()
ax[0].set_xlabel(r"$\varepsilon$",fontsize = 24)
ax[0].set_ylabel(r"${\sigma _{tt}}$, ГПа",fontsize = 24)
ax[0].legend(loc='best')

points4= np.genfromtxt("snn_a100_Ms60_1.csv", delimiter=",")
points5= np.genfromtxt("snn_a100_Ms60_T0.csv", delimiter=",")
points6= np.genfromtxt("snn_a100_Ms60_sig0.csv", delimiter=",")

ax[1].plot(points4[0, 4:10], points4[-1, 4:10] /1e9,'r',label="$T=0.1$ ГПа; ${\sigma _{0}^{s}}=1$ Н/м")
ax[1].plot(points5[0, 4:10], points5[-1, 4:10] /1e9,'g',label="$T=0$ ГПа; ${\sigma _{0}^{s}}=1$ Н/м" )
ax[1].plot(points6[0, 4:10], points6[-1, 4:10] /1e9, 'b',label="$T=0.1$ ГПа; ${\sigma _{0}^{s}}=0$ Н/м")

ax[1].plot(points4[0, 4:10], points4[1, 4:10] /1e9,'r--')
ax[1].plot(points5[0, 4:10], points5[1, 4:10] /1e9,'g--')
ax[1].plot(points6[0, 4:10], points6[1, 4:10] /1e9,'b--')

#ax[1].annotate(r"$1$", xy=(13.5,0.052), xytext=(15,0.04),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)
#ax[1].annotate(r"$1$", xy=(14,0.0333), xytext=(15,0.04),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)

#ax[1].annotate(r"$2$", xy=(14,0.0557), xytext=(15,0.054),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)
#ax[1].annotate(r"$2$", xy=(14,0.0444), xytext=(15,0.054),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)

#ax[1].annotate(r"$3$", xy=(14,0), xytext=(15,-0.005),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)
#ax[1].annotate(r"$3$", xy=(14,-0.0109), xytext=(15,-0.005),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)

#leg = ax[1].legend(['1 - $T=0.1$ GPa; ${\sigma _{0}^{s}}=1$ N/m','2 - $T=0$ GPa; ${\sigma _{0}^{s}}=1$ N/m', '3 - $T=0.1$ GPa; ${\sigma _{0}^{s}}=0$ N/m'],loc='best',handlelength=0, handletextpad=0, fancybox=True)
#for item in leg.legendHandles:
 #   item.set_visible(False)

ax[1].set_xticks(points1[0,4:10:2], [ '%.2f' % elem for elem in points1[0,4:10:2] ])

ax[1].grid()
ax[1].set_xlabel(r"$\varepsilon$",fontsize = 24)
ax[1].set_ylabel(r"${\sigma _{nn}}$, ГПа",fontsize = 24)
ax[1].legend(loc='best')

plt.show()