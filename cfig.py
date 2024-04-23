import numpy as np
from matplotlib import pyplot as plt

points1= np.genfromtxt("stt_n15.csv", delimiter=",")
points2= np.genfromtxt("stt_n15_T0.csv", delimiter=",")
points3= np.genfromtxt("stt_n15_sig0.csv", delimiter=",")

a  = 10

fig, ax = plt.subplots(nrows=1, ncols=2,layout="constrained")
ax[0].tick_params(axis='both', which='major', labelsize=18)
ax[1].tick_params(axis='both', which='major', labelsize=18)

ax[0].plot(points1[0, 4:15]*a, points1[-1, 4:15] /1e9,'r',label="$T=0.1$ GPa; ${\sigma _{0}^{s}}=1$ N/m")
ax[0].plot(points2[0, 4:15]*a, points2[-1, 4:15] /1e9,'g',label="$T=0$ GPa; ${\sigma _{0}^{s}}=1$ N/m" )
ax[0].plot(points3[0, 4:15]*a, points3[-1, 4:15] /1e9, 'b',label="$T=0.1$ GPa; ${\sigma _{0}^{s}}=0$ N/m")

ax[0].plot(points1[0, 4:15]*a, points1[1, 4:15] /1e9,'r--')
ax[0].plot(points2[0, 4:15]*a, points2[1, 4:15] /1e9,'g--')
ax[0].plot(points3[0, 4:15]*a, points3[1, 4:15] /1e9,'b--')


#ax[0].annotate(r"$1$", xy=(1.27,0.042), xytext=(1.33,0.37),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)
#ax[0].annotate(r"$1$", xy=(1.21,0.69), xytext=(1.33,0.37),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)

#ax[0].annotate(r"$2$", xy=(1.34,-0.14), xytext=(1.49,0.43),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)
#ax[0].annotate(r"$2$", xy=(1.44,0.53), xytext=(1.49,0.43),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)

#ax[0].annotate(r"$3$", xy=(1.5,0.15), xytext=(1.52,0.21),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)
#ax[0].annotate(r"$3$", xy=(1.5,0.27), xytext=(1.52,0.21),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)

#leg = ax[0].legend(['1 - $T=0.1$ GPa; ${\sigma _{0}^{s}}=1$ N/m','2 - $T=0$ GPa; ${\sigma _{0}^{s}}=1$ N/m', '3 - $T=0.1$ GPa; ${\sigma _{0}^{s}}=0$ N/m'],loc='best',handlelength=0, handletextpad=0, fancybox=True)
#for item in leg.legendHandles:
#    item.set_visible(False)
#plt.title("snn a=100")


#xi = list(range(len(points1[0,4:15]*a)))
ax[0].set_xticks(points1[0,4:15:2]*a, [ '%.1f' % elem for elem in points1[0,4:15:2]*a ])


ax[0].grid()
ax[0].set_xlabel(r"$A$, nm",fontsize = 24)
ax[0].set_ylabel(r"${\sigma _{tt}}$, GPa",fontsize = 24)
ax[0].legend(loc='best')


points4= np.genfromtxt("snn_n15.csv", delimiter=",")
points5= np.genfromtxt("snn_n15_T0.csv", delimiter=",")
points6= np.genfromtxt("snn_n15_sig0.csv", delimiter=",")

ax[1].plot(points4[0, 4:15]*a, points4[-1, 4:15] /1e9,'r',label="$T=0.1$ GPa; ${\sigma _{0}^{s}}=1$ N/m")
ax[1].plot(points5[0, 4:15]*a, points5[-1, 4:15] /1e9,'g',label="$T=0$ GPa; ${\sigma _{0}^{s}}=1$ N/m" )
ax[1].plot(points6[0, 4:15]*a, points6[-1, 4:15] /1e9, 'b',label="$T=0.1$ GPa; ${\sigma _{0}^{s}}=0$ N/m")

ax[1].plot(points4[0, 4:15]*a, points4[1, 4:15] /1e9,'r--')
ax[1].plot(points5[0, 4:15]*a, points5[1, 4:15] /1e9,'g--')
ax[1].plot(points6[0, 4:15]*a, points6[1, 4:15] /1e9,'b--')

#ax[1].annotate(r"$1$", xy=(1.243,0.584), xytext=(1.142,0.587),arrowprops=dict(arrowstyle="-",facecolor='black'),fontsize=20)
#ax[1].annotate(r"$1$", xy=(1.233,0.486), xytext=(1.142,0.587),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)

#ax[1].annotate(r"$2$", xy=(1.437,0.636), xytext=(1.5,0.618),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)
#ax[1].annotate(r"$2$", xy=(1.423,0.555), xytext=(1.5,0.618),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)

#ax[1].annotate(r"$3$", xy=(1.4,0.032), xytext=(1.5,0.06),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)
#ax[1].annotate(r"$3$", xy=(1.4,0.005), xytext=(1.5,0.06),arrowprops=dict(arrowstyle="-",facecolor='black'), fontsize=20)

#leg = ax[1].legend(['1 - $T=0.1$ GPa; ${\sigma _{0}^{s}}=1$ N/m','2 - $T=0$ GPa; ${\sigma _{0}^{s}}=1$ N/m', '3 - $T=0.1$ GPa; ${\sigma _{0}^{s}}=0$ N/m'],loc='best',handlelength=0, handletextpad=0, fancybox=True)
#for item in leg.legendHandles:
#    item.set_visible(False)

ax[1].set_xticks(points1[0,4:15:2]*a, [ '%.1f' % elem for elem in points1[0,4:15:2]*a ])

ax[1].grid()
ax[1].set_xlabel(r"$A$, nm",fontsize = 24)
ax[1].set_ylabel(r"${\sigma _{nn}}$, GPa",fontsize = 24)
ax[1].legend(loc='best')

plt.show()