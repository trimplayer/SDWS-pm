import numpy as np
from matplotlib import pyplot as plt

points1= np.genfromtxt("stt_a10_Ms6_1.csv", delimiter=",")

for i in range(1,11):
    plt.plot(points1[0, :], points1[i, :]/1e9, label='$n={n}$'.format(n=i))


for j in range(1,10):
    p = 0
    n1 = j
    n2 = j+1
    print(points1[0, p])
    #print(points1[n1, p])
    #print(points1[n2, p])
    print(n2)
    print((points1[n2, p]-points1[n1, p])/points1[n1, p])
    print((points1[n2, p] - points1[1, p]) / points1[1, p])


plt.legend(loc='best')
plt.grid()
plt.xlabel(r"$\varepsilon$",fontsize = 20)
plt.ylabel(r"${\sigma _{nn}}$, ГПа",fontsize = 20)
plt.show()