import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt

# whole
in_file = open("points1","rb")
points1 = pickle.load(in_file)
in_file.close()

in_file = open("points2","rb")
points2 = pickle.load(in_file)
in_file.close()

in_file = open("points3","rb")
points3 = pickle.load(in_file)
in_file.close()

in_file = open("points4","rb")
points4 = pickle.load(in_file)
in_file.close()

in_file = open("points5","rb")
points5 = pickle.load(in_file)
in_file.close()

in_file = open("points6","rb")
points6 = pickle.load(in_file)
in_file.close()

in_file = open("points7","rb")
points7 = pickle.load(in_file)
in_file.close()

in_file = open("points8","rb")
points8 = pickle.load(in_file)
in_file.close()

fig, ax = plt.subplots(nrows=1, ncols=2)

ax[0].plot(points1[:, 0], points1[:, 1], label="e=0.1", color="black")
ax[0].plot(points2[:, 0], points2[:, 1], color="black")
ax[0].plot(points5[:, 0], points5[:, 1],"--", label="e=0.05", color="black")
ax[0].plot(points6[:, 0], points6[:, 1],"--", color="black")

ax[0].annotate(r"$n = 1$", xy=(0,0.33), xytext=(0.3,0.3),arrowprops=dict(arrowstyle="-",facecolor='black'))
ax[0].annotate(r"$n = 1$", xy=(0,0.58), xytext=(0.3,0.3),arrowprops=dict(arrowstyle="-",facecolor='black'))

ax[0].annotate(r"$n = 2$", xy=(0,-0.05), xytext=(-0.4,0.3),arrowprops=dict(arrowstyle="-",facecolor='black'))
ax[0].annotate(r"$n = 2$", xy=(0,0.2), xytext=(-0.4,0.3),arrowprops=dict(arrowstyle="-",facecolor='black'))

ax[0].set_xlim([-0.49,0.49])
ax[0].set_ylim([-0.7,0.7])

ax[0].grid()
ax[0].set_xlabel(r"${x_1}/a$")
ax[0].set_ylabel(r"${\sigma _{tt}}$, ГПа")
ax[0].legend( prop={'size': 12})

ax[1].plot(points3[:, 0], points3[:, 1], label="e=0.1", color="black")
ax[1].plot(points4[:, 0], points4[:, 1], color="black")
ax[1].plot(points7[:, 0], points7[:, 1],"--", label="e=0.05", color="black")
ax[1].plot(points8[:, 0], points8[:, 1],"--", color="black")

ax[1].annotate(r"$n = 1$", xy=(0,0.38), xytext=(0.3,0.2),arrowprops=dict(arrowstyle="-",facecolor='black'))
ax[1].annotate(r"$n = 1$", xy=(0,0.18), xytext=(0.3,0.2),arrowprops=dict(arrowstyle="-",facecolor='black'))

ax[1].annotate(r"$n = 2$", xy=(0,0.52), xytext=(-0.3,0.4),arrowprops=dict(arrowstyle="-",facecolor='black'))
ax[1].annotate(r"$n = 2$", xy=(0,0.25), xytext=(-0.3,0.4),arrowprops=dict(arrowstyle="-",facecolor='black'))

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

