from values import *
from sympy import *
from func import ser

n = 1
e = symbols('e')
f = symbols('f')


phip1i1 = IndexedBase('Phip1_i_1')
phip1ij1 = IndexedBase('Phi_p1_ij_1')
i,j = symbols('i j', cls=Idx)

phip11 = Sum(e**i/factorial(i) * phip1i1[i], (i,0,n)).doit()

for i in range(n+1):
    phip11 = phip11.subs(phip1i1[i], Sum((I*e*f)**j / factorial(j) * phip1ij1[i, j],(j,0,n)).doit())

phip11 = ser(phip11, n+1)

out_file = open("phi1p1.txt","w")
out_file.write(str(phip11))
out_file.close()

dphip1i1 = IndexedBase('dPhip1_i_1')
i,j = symbols('i j', cls=Idx)
dphip11 = Sum(e**i/factorial(i) * dphip1i1[i], (i,0,n)).doit()

for i in range(n+1):
    dphip11 = dphip11.subs(dphip1i1[i], Sum((I*e*f)**j / factorial(j) * phip1ij1[i, j+1],(j,0,n)).doit())

dphip11 = ser(dphip11, n+1)

out_file = open("dphi1p1.txt","w")
out_file.write(str(dphip11))
out_file.close()

Upsilonp1i1 = IndexedBase('Upsilonp1_i_1')
Upsilonp1ij1 = IndexedBase('Upsilon_p1_ij_1')
i,j = symbols('i j', cls=Idx)

Upsilonp11 = Sum(e**i/factorial(i) * Upsilonp1i1[i], (i,0,n)).doit()

for i in range(n+1):
    Upsilonp11 = Upsilonp11.subs(Upsilonp1i1[i], Sum((I*e*f)**j / factorial(j) * Upsilonp1ij1[i, j],(j,0,n)).doit())

Upsilonp11 = ser(Upsilonp11, n+1)

out_file = open("Upsilon1p1.txt","w")
out_file.write(str(Upsilonp11))
out_file.close()

