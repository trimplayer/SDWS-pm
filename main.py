from values import *
from sympy import *
from sympy.polys.ring_series import rs_series
from func import ser

n = 1
e = symbols('e')
f = symbols('f')
df = symbols('df')
d2f = symbols('d2f')

duu = []

for i in range(n+1):
    duu.append(symbols('dus' + str(i)))

dus = Float(0)
for i in range(n+1):
    dus += (e**i/factorial(i)) * duu[i]

duuu = []
for i in range(n+1):
    duu1 = []
    for j in range(n+2):
        duu1.append(symbols('dus' + str(i) + str(j)))
    duuu.append(duu1)


for i in range(n+1):
    sub_ser = Float(0)
    for j in range(n+1):
        sub_ser += ((I*e*f)**j/factorial(j)) * duuu[i][j]
    dus = dus.subs(duu[i], sub_ser)

dus = ser(dus, n+1)

d2uu = []

for i in range(n+1):
    d2uu.append(symbols('dus' + str(i)))\

d2us = Float(0)
for i in range(n+1):
    d2us += (e**i/factorial(i)) * d2uu[i]

for i in range(n+1):
    sub_ser = Float(0)
    for j in range(n+1):
        sub_ser += ((I*e*f)**j/factorial(j)) * duuu[i][j+1]
    d2us = d2us.subs(d2uu[i], sub_ser)

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
