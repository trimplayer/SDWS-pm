import numpy as np
from matplotlib import pyplot as plt
import csv

points1= np.genfromtxt("snn_a100_Ms60_1.csv", delimiter=",")
points2= np.genfromtxt("snn_a100_Ms60_T0.csv", delimiter=",")
points3= np.genfromtxt("snn_a100_Ms60_sig0.csv", delimiter=",")
points4= np.genfromtxt("snn_a100_Ms60_psms.csv", delimiter=",")
points5= np.genfromtxt("snn_a100_Ms60_psms_T0.csv", delimiter=",")
filename = "rc_snn_a100_Ms60.csv"

casedic = {id(points1):'1',id(points2):'T = 0',id(points3):'sig = 0',id(points4):'ps = ms',id(points5):'ps = ms T = 0'}

fields = ['n','case','eps', 'rc to prev', 'rc to 1']

# for j in range(1,10):
#     p = 4
#     n1 = j
#     n2 = j+1
#     print("eps = ",points1[0, p]," n1 = ",n1," n2 = ",n2)
#     #print(points[n1, p])
#     #print(points[n2, p])
#     print((points1[n2, p]-points1[n1, p])/points1[n1, p])

plist = [points1, points2, points3, points4, points5]
#epslist = [0.05,0.1,0.15,0.2]
epslist = [4,9,14]

frows =  []

for pi in plist:
    for epsi in epslist:
        p = epsi
        for j in range(1, 10):
            n1 = j
            n2 = j + 1
            #print("eps = ", pi[0, p], " n1 = ", n1, " n2 = ", n2)
            #print((pi[n2, p] - pi[n1, p]) / pi[n1, p])
            rel_cha = (pi[n2, p] - pi[n1, p]) / pi[n1, p]
            rel_cha_to1 = (pi[n2, p] - pi[1, p]) / pi[1, p]
            row1 = [n2,casedic[id(pi)],points1[0, p],rel_cha,rel_cha_to1]
            frows.append(row1)
            print(row1)

with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)

    # writing the fields
    csvwriter.writerow(fields)

    # writing the data rows
    csvwriter.writerows(frows)
