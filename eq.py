from values import *
from sympy import *
from func import ser

n = 1
e = symbols('e')
f = symbols('f')


phip1i1 = IndexedBase('Phip1_i_1')
phip1ij1 = IndexedBase('Phip1_ij_1')
i,j = symbols('i j', cls=Idx)

phip11 = Sum(e**i/factorial(i) * phip1i1[i], (i,0,n)).doit()

for i in range(n+1):
    phip11 = phip11.subs(phip1i1[i], Sum((I*e*f)**j / factorial(j) * phip1ij1[i, j],(j,0,n)).doit())

phip11 = ser(phip11, n+1)
