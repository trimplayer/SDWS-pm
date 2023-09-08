import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt


# T = 0
in_file = open("points17","rb")
points17 = pickle.load(in_file)
in_file.close()

in_file = open("points18","rb")
points18 = pickle.load(in_file)
in_file.close()

in_file = open("points19","rb")
points19 = pickle.load(in_file)
in_file.close()

in_file = open("points20","rb")
points20 = pickle.load(in_file)
in_file.close()

in_file = open("points21","rb")
points21 = pickle.load(in_file)
in_file.close()

in_file = open("points22","rb")
points22 = pickle.load(in_file)
in_file.close()

in_file = open("points23","rb")
points23 = pickle.load(in_file)
in_file.close()

in_file = open("points24","rb")
points24 = pickle.load(in_file)
in_file.close()

fig, ax = plt.subplots(nrows=1, ncols=2)

ax[0].plot(points17[:, 0], points17[:, 1], label="n=1 e=0.1")
ax[0].plot(points18[:, 0], points18[:, 1], label="n=2 e=0.1")
ax[0].plot(points21[:, 0], points21[:, 1], label="n=1 e=0.05")
ax[0].plot(points22[:, 0], points22[:, 1], label="n=2 e=0.05")

ax[0].set_xlim([-0.49,0.49])
ax[0].set_ylim([-0.7,0.7])

ax[0].grid()
ax[0].set_xlabel(r"${x_1}/a$")
ax[0].set_ylabel(r"${\sigma _{tt}}$, ГПа")
ax[0].legend( prop={'size': 12})

ax[1].plot(points19[:, 0], points19[:, 1], label="n=1 e=0.1")
ax[1].plot(points20[:, 0], points20[:, 1], label="n=2 e=0.1")
ax[1].plot(points23[:, 0], points23[:, 1], label="n=1 e=0.05")
ax[1].plot(points24[:, 0], points24[:, 1], label="n=2 e=0.05")

ax[1].set_xlim([-0.49,0.49])
ax[1].set_ylim([-0.7,0.7])

ax[1].grid()
ax[1].set_xlabel(r"${x_1}/a$")
ax[1].set_ylabel(r"${\sigma _{nn}}$, ГПа")
ax[1].legend(prop={'size': 12})

for item in ([ax[0].title, ax[0].xaxis.label, ax[0].yaxis.label]):
    item.set_fontsize(18)
for item in ([ax[1].title, ax[1].xaxis.label, ax[1].yaxis.label]):
    item.set_fontsize(18)

plt.show()

#plt.savefig('fig02.pdf', format='pdf', dpi=300, bbox_inches="tight")
#plt.savefig('fig02.jpg', format='jpg', dpi=300, bbox_inches="tight")

