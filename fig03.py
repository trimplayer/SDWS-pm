import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt


# sig = 0
in_file = open("points9","rb")
points9 = pickle.load(in_file)
in_file.close()

in_file = open("points10","rb")
points10 = pickle.load(in_file)
in_file.close()

in_file = open("points11","rb")
points11 = pickle.load(in_file)
in_file.close()

in_file = open("points12","rb")
points12 = pickle.load(in_file)
in_file.close()

in_file = open("points13","rb")
points13 = pickle.load(in_file)
in_file.close()

in_file = open("points14","rb")
points14 = pickle.load(in_file)
in_file.close()

in_file = open("points15","rb")
points15 = pickle.load(in_file)
in_file.close()

in_file = open("points16","rb")
points16 = pickle.load(in_file)
in_file.close()

fig, ax = plt.subplots(nrows=1, ncols=2)

ax[0].plot(points9[:, 0], points9[:, 1], label="n=1 e=0.1")
ax[0].plot(points10[:, 0], points10[:, 1], label="n=2 e=0.1")
ax[0].plot(points13[:, 0], points13[:, 1], label="n=1 e=0.05")
ax[0].plot(points14[:, 0], points14[:, 1], label="n=2 e=0.05")

ax[0].set_xlim([-0.49,0.49])
ax[0].set_ylim([-0.7,0.7])

ax[0].grid()
ax[0].set_xlabel(r"${x_1}/a$")
ax[0].set_ylabel(r"${\sigma _{tt}}$, ГПа")
ax[0].legend( prop={'size': 12})

ax[1].plot(points11[:, 0], points11[:, 1], label="n=1 e=0.1")
ax[1].plot(points12[:, 0], points12[:, 1], label="n=2 e=0.1")
ax[1].plot(points15[:, 0], points15[:, 1], label="n=1 e=0.05")
ax[1].plot(points16[:, 0], points16[:, 1], label="n=2 e=0.05")

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

