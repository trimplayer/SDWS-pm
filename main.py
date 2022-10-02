from values import *
from sympy import *

n = 1
e = symbols('e')
f = symbols('f')

duu = []

for i in range(n+1):
    duu.append(symbols('dus' + str(i)))

dus = Float(0)
for i in range(n+1):
    dus += (e**i/factorial(i)) * duu[i]
#dus -= e

#pprint(dus)

#(I*e*f)^j/factorial(j)

duuu = []
for i in range(n+1):
    duu1 = []
    for j in range(n+1):
        duu1.append(symbols('dus' + str(i) + str(j)))
    duuu.append(duu1)


for i in range(n+1):
    sub_ser = Float(0)
    for j in range(n+1):
        sub_ser += ((I*e*f)**j/factorial(j)) * duuu[i][j]
    dus = dus.subs(duu[i], sub_ser)

