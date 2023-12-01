import numpy as np
from matplotlib import pyplot as plt

points1= np.genfromtxt("stt_n10.csv", delimiter=",")

for i in range(1,11):
    plt.plot(points1[0, :], points1[i, :], label='$n={n}$'.format(n=i))

plt.legend(loc='best')
plt.show()