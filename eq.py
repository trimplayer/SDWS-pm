from values import *
from sympy import *
from func import *
import pickle

n = 2
e = symbols('e', real=True)
f = symbols('f', real=True)
df = symbols('df', real=True)
d2f = symbols('d2f', real=True)
d = {'e' : e, 'f' : f, 'df' : df, 'd2f' : d2f}


phip1i1 = IndexedBase('Phip1_i_1')
phip1ij1 = IndexedBase('Phi_p1_1')
i,j = symbols('i j', cls=Idx)

phip11 = Sum(e**i/factorial(i) * phip1i1[i], (i,0,n)).doit()

for i in range(n+1):
    phip11 = phip11.subs(phip1i1[i], Sum((I*e*f)**j / factorial(j) * phip1ij1[i, j],(j,0,n)).doit())

phip11 = ser(phip11, n+1)

out_file = open("phi1p1","wb")
pickle.dump(phip11,out_file)
out_file.close()

dphip1i1 = IndexedBase('dPhip1_i_1')
i,j = symbols('i j', cls=Idx)
dphip11 = Sum(e**i/factorial(i) * dphip1i1[i], (i,0,n)).doit()

for i in range(n+1):
    dphip11 = dphip11.subs(dphip1i1[i], Sum((I*e*f)**j / factorial(j) * phip1ij1[i, j+1],(j,0,n)).doit())

dphip11 = ser(dphip11, n+1)

out_file = open("dphi1p1","wb")
pickle.dump(dphip11, out_file)
out_file.close()

Upsilonp1i1 = IndexedBase('Upsilonp1_i_1')
Upsilonp1ij1 = IndexedBase('Upsilon_p1_1')
i,j = symbols('i j', cls=Idx)

Upsilonp11 = Sum(e**i/factorial(i) * Upsilonp1i1[i], (i,0,n)).doit()

for i in range(n+1):
    Upsilonp11 = Upsilonp11.subs(Upsilonp1i1[i], Sum((-I*e*f)**j / factorial(j) * Upsilonp1ij1[i, j],(j,0,n)).doit())

Upsilonp11 = ser(Upsilonp11, n+1)

out_file = open("Upsilon1p1","wb")
pickle.dump(Upsilonp11, out_file)
out_file.close()

expoa = 1 - 2*I*e*df/(1+I*e*df)
expoa = expoa.series(e, n=n+1)
expoa = expoa.removeO()

in_file = open("rhs1","rb")
rhs1 = pickle.load(in_file)
in_file.close()

w1 = symbols('w1')
x = symbols('x', real=True)

eqq1 = phip11 + conjugate(phip11) - (Upsilonp11 + conjugate(phip11) - (w1 - conjugate(w1)) * conjugate(dphip11)) * expoa - rhs1
eqq1 = eqq1.subs(w1,x+I*e*f)

eqq1 = ser(eqq1, n+1)



eqq1l =[]
for l in range(n+1):
    eqq1l.append(eqq1.coeff(e, l))

#eqq11 = eqq1l[1]  ##
eqq11 = eqq1

out_file = open("eqq11","wb")
pickle.dump(eqq11, out_file)
out_file.close()
