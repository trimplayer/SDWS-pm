from values import *
from sympy import *
from func import *
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pickle

n = 1

f = symbols('f', real=True)
df = symbols('df', real=True)
d2f = symbols('d2f', real=True)
e = symbols('e', real=True)
ka = symbols('kappa', real=True)
mu = symbols('mu', real=True)
x = symbols('x', real=True)
A = IndexedBase('A', real=True)
B = IndexedBase('B', real=True)

unkn=[A[1, 1],A[1, (-1)],B[1, 1],B[1, (-1)]]


in_file = open("cfc12","rb")
cfc12 = pickle.load(in_file)
in_file.close()

in_file = open("cfs12","rb")
cfs12 = pickle.load(in_file)
in_file.close()

eqs = [re(cfc12), im(cfc12), re(cfs12), im(cfs12)]

qs0 = symbols('q_s0', real=True) ##
par0 = {qs0: 0}

nu = symbols("nu")
lam = symbols("lambda")
par1 = {b: 2*pi /a, ka: 3-4*(lam/(2*(lam+mu))), x: x*a, nu: lam/(2*(lam+mu))}

sigmass0 = symbols('sigma_ss0')
Ms = symbols("M_s")
sigmas0 = symbols("sigma_s0")
T = symbols('T', real=True)
par2 = {lam:58.17*10**9, mu:26.13*10**9, Ms:6.099, a:10*10**(-9), T:0.1*10**9,sigmass0:0,sigmas0:0}

i = 0
for ex in eqs:
    eqs[i] = ex.subs(par0)
    eqs[i] = eqs[i].subs(par1)
    eqs[i] = eqs[i].evalf(subs=par2)
    i +=1

sol = solve(eqs, unkn)
subc = sol

in_file = open("Upsilon1p1N","rb")
Upsilon1p1 = pickle.load(in_file)
in_file.close()

in_file = open("phi1p1N","rb")
phi1p1 = pickle.load(in_file)
in_file.close()

in_file = open("dphi1p1N","rb")
dphi1p1 = pickle.load(in_file)
in_file.close()

subf = {f:sub_f, df:sub_df, d2f:sub_d2f}

z = symbols('z')

Upsilon1p1 = Upsilon1p1.subs(z, x)
Upsilon1p1 = Upsilon1p1.subs(subc)
Upsilon1p1 = Upsilon1p1.subs(par0)
Upsilon1p1 = Upsilon1p1.subs(par1)
Upsilon1p1 = Upsilon1p1.subs(par2)

phi1p1 = phi1p1.subs(z, x)
phi1p1 = phi1p1.subs(subc)
phi1p1 = phi1p1.subs(par0)
phi1p1 = phi1p1.subs(par1)
phi1p1 = phi1p1.subs(par2)
phi1p1 = phi1p1.evalf()

dphi1p1 = dphi1p1.subs(z, x)
dphi1p1 = dphi1p1.subs(subc)
dphi1p1 = dphi1p1.subs(par0)
dphi1p1 = dphi1p1.subs(par1)
dphi1p1 = dphi1p1.subs(par2)
dphi1p1 = dphi1p1.evalf()

expoa = 1 - 2 * I * e * df / (1 + I * e * df)
expoa = expoa.series(e, n=n+1)
expoa = expoa.removeO()

h = (1 + e**2 * df**2)
hm = series(1/h ,e ,n=n+1)
hm = hm.removeO()

cr = e * d2f / h ** 3
cr = cr.series(e, n=n+1)
cr = cr.removeO()

w1 = symbols("w1", complex=True)
G1 = phi1p1+conjugate(phi1p1)-(Upsilon1p1+conjugate(phi1p1)-(w1-conjugate(w1))*conjugate(dphi1p1))*expoa
G1 = G1.subs(w1, x+I*e*f)

G1 = ser(G1, n+1)

sigma_1nn = re(G1)

sigma_1tt_plus_sigma_1nn = re(4*phi1p1)


sigma_1tt=sigma_1tt_plus_sigma_1nn-sigma_1nn  ##

sigmann = sigma_1nn
sigmann = sigmann.subs(e, 0.1)
scf = sigmann.evalf(subs={x: 0})

plot(sigmann, (x, -0.5, 0.5))

