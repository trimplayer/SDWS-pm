from values import *
from sympy import *
from func import *
import pickle

n = 4
e = symbols('e', real=True)
f = symbols('f', real=True)
df = symbols('df', real=True)
d2f = symbols('d2f', real=True)

dusi = IndexedBase("dus_i")
dusij = IndexedBase("du_s")
i, j = symbols('i j', cls=Idx)

dus = Sum((e ** i / factorial(i)) * dusi[i], (i, 0, n)).doit()

for k in range(n + 1):
    dus = dus.subs(dusi[k], Sum(((I * e * f) ** j / factorial(j)) * dusij[k, j], (j, 0, n)).doit())

dus = ser(dus,n+1)


d2usi = IndexedBase("dus_i")

d2us = Sum((e**i/factorial(i)) * d2usi[i],(i,0,n)).doit()

for k in range(n + 1):
    d2us = d2us.subs(d2usi[k],Sum(((I*e*f)**j/factorial(j)) * dusij[k, j+1], (j,0,n)).doit())

d2us = ser(d2us, n+1)

expoa = (1 + I*e*df) / sqrt(1 + e**2 * df**2)
expoa = expoa.series(e, n=n+1)
expoa = expoa.removeO()


h = (1 + e**2 * df**2)  #
hm = series((I + e * df)/h ,e ,n=n+1)
hm = hm.removeO()

cr = e * d2f / h ** 3
cr = cr.series(e, n=n+1)
cr = cr.removeO()

sigmass0 = symbols('sigma_ss0')
Ms = symbols("M_s")
sigmas0 = symbols("sigma_s0")

rhs1 = (sigmass0 * cr + Ms * cr * re(dus) + sigmas0 * im(d2us * expoa)) - I * (Ms * re(d2us * expoa) - sigmas0 * cr * im(dus))
rhs1 = ser(rhs1, n+1)

print(rhs1)

out_file = open("rhs1","wb")
pickle.dump(rhs1,out_file)
out_file.close()