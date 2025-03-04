import numpy as np
from matplotlib import pyplot as plt

points1= np.genfromtxt("stt_n15_sig0.csv", delimiter=",")
points2 = np.genfromtxt("sigtt_nlg_conhe_sig0.csv")

points3= np.genfromtxt("stt_n15.csv", delimiter=",")
points4 = np.genfromtxt("sigtt_nlg_conhe.csv")
points5=np.genfromtxt("sig.csv")
#points3=np.genfromtxt("sigtt_conhe.csv")
#points4=np.genfromtxt("sigtt_nlg.csv")
#points5=np.genfromtxt("sigtt_nlg_conhe.csv")




#for i in range(1,11):
#    plt.plot(points1[0, :], points1[i, :]/1e9, label='$n={n}$'.format(n=i))

plt.plot(points1[0,:],points1[1,:]/1e9,label='n=1 sig=0')
plt.plot(points1[0,:],points1[15,:]/1e9,label='n=15 sig =0')
plt.plot(points2[:,0], points2[:,1]*1e9,label='fe sig=0')

plt.plot(points3[0,:],points3[1,:]/1e9,label='n=1')
plt.plot(points3[0,:],points3[15,:]/1e9,label='n=15')
plt.plot(points4[:,0], points4[:,1]*1e9,label='fe')

plt.plot(points5[:,0], points5[:,1]*1e9,label='fe init')

#plt.plot(points2[:,0], points2[:,1]*1e9,label='fe')
#plt.plot(points3[:,0], points3[:,1]*1e9,label='fe constant height')
#plt.plot(points4[:,0], points4[:,1]*1e9,label='fe nlgeom on')
#plt.plot(points5[:,0], points5[:,1]*1e9,label='fe nlgeom on con h')

#plt.ylim([-0.1,0.4])

plt.plot()

for j in range(1,15):
    p = 14
    n1 = j
    n2 = j+1
    print(points1[0, p])
    #print(points1[n1, p])
    #print(points1[n2, p])
    print(n2)
    print(abs(points1[n2, p]-points1[n1, p])/abs(points1[n1, p]))
    print(abs(points1[n2, p] - points1[1, p]) / abs(points1[1, p]))


plt.legend(loc='best')
plt.grid()
plt.xlabel(r"$\varepsilon$",fontsize = 20)
plt.ylabel(r"${\sigma _{tt}}$, ГПа",fontsize = 20)
plt.title("M=6 sig_s=0")
plt.tight_layout()
plt.show()