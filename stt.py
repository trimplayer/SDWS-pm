from values import *
from sympy import *
from func import *
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pickle

n = 4

f = symbols('f', real=True)
df = symbols('df', real=True)
d2f = symbols('d2f', real=True)
e = symbols('e', real=True)
ka = symbols('kappa', real=True)
mu = symbols('mu', real=True)
x = symbols('x', real=True)
A = IndexedBase('A', real=True)
B = IndexedBase('B', real=True)

#unkn=[A[1, 1],A[1, (-1)],B[1, 1],B[1, (-1)]]
#unkn2=[A[2, 1],A[2, (-1)],B[2, 1],B[2, (-1)]]
unkn = [0]

for k in range(1,n+1):
    unkn.append([A[k, 1],A[k, (-1)],B[k, 1],B[k, (-1)]])

# in_file = open("cfc12","rb")
# cfc12 = pickle.load(in_file)
# in_file.close()
#
# in_file = open("cfs12","rb")
# cfs12 = pickle.load(in_file)
# in_file.close()

# eqs = [re(cfc12), im(cfc12), re(cfs12), im(cfs12)]

in_file = open("eqslist","rb")
eqslist = pickle.load(in_file)
in_file.close()



qs0 = symbols('q_s0', real=True) ##
par0 = {qs0: 0}

nu = symbols("nu")
lam = symbols("lambda")
par1 = {b: 2*pi /a, ka: 3-4*(lam/(2*(lam+mu))), x: x*a, nu: lam/(2*(lam+mu))}

sigmass0 = symbols('sigma_ss0')
Ms = symbols("M_s")
sigmas0 = symbols("sigma_s0")
T = symbols('T', real=True)
par2 = {lam:58.17*10**9, mu:26.13*10**9, Ms:6.099, a:10*10**(-9), T:0.1*10**9,sigmass0:1,sigmas0:1} #origin
#par2 = {lam:58.17*10**9, mu:26.13*10**9, Ms:6.099, a:10*10**(-9), T:0.1*10**9,sigmass0:0,sigmas0:0}

# eq2 = -1/16/mu*re(-8*((-a*sigmass0*b+(I*B[1, 1]+A[1, 1])*Ms-T*a+sigmas0*(I*A[1, -1]-B[1, -1]))*mu-1/8*a*b*T*Ms*(ka+1))*(ka+1)*b-(I*32*A[1, -1]+32*A[1, 1])*mu**2-16*a*b**2*sigmass0*mu-2*a*b**2*Ms*T*(ka+1))
# eq4 = -1/16/mu*im(-8*((-a*sigmass0*b+(I*B[1, 1]+A[1, 1])*Ms-T*a+sigmas0*(I*A[1, -1]-B[1, -1]))*mu-1/8*a*b*T*Ms*(ka+1))*(ka+1)*b-(I*32*A[1, -1]+32*A[1, 1])*mu**2-16*a*b**2*sigmass0*mu-2*a*b**2*Ms*T*(ka+1))
# eq6 = 1/2/mu*re(-I*((-a*sigmass0*b+(I*B[1, 1]+A[1, 1])*Ms-T*a+sigmas0*(I*A[1, -1]-B[1, -1]))*mu-1/8*a*b*T*Ms*(ka+1))*(ka+1)*b+4*mu**2*(B[1, 1]+I*B[1, -1]))
# eq8 = 1/2/mu*im(-I*((-a*sigmass0*b+(I*B[1, 1]+A[1, 1])*Ms-T*a+sigmas0*(I*A[1, -1]-B[1, -1]))*mu-1/8*a*b*T*Ms*(ka+1))*(ka+1)*b+4*mu**2*(B[1, 1]+I*B[1, -1]))
#
# eqs = [eq2, eq4, eq6, eq8]

#i = 0
# for ex in eqs:
#     eqs[i] = ex.subs(par0)
#     eqs[i] = eqs[i].subs(par1)
#     #eqs[i] = eqs[i].subs(par1)
#     eqs[i] = eqs[i].evalf(subs=par2)
#     #eqs[i] = eqs[i].subs(par2)
#     i +=1
#
#
# sol = solve(eqs, unkn)
# subc = sol

# unkn=[A[1, 1],A[1, (-1)],B[1, 1],B[1, (-1)]]
# unkn2=[A[2, 1],A[2, (-1)],B[2, 1],B[2, (-1)]]
# eqss = eqslist[0].copy()
# i = 0
# for ex in eqss:
#     eqss[i] = ex.subs(par0)
#     eqss[ i] = eqss[ i].subs(par1)
#     eqss[ i] = eqss[ i].evalf(subs=par2)
#     i += 1
# sol1 = solve(eqss, unkn)
# subc1 = sol1
#
# i = 0
# eqs0 = [0,0,0,0]
# for ex in eqslist[3]:
#     eqs0[i]=ex.copy()
#     # eqs0[i] = eqs0[i] + eqslist[2][i]  # eqs0[2,1] - sin(2*pix*x)
#     eqs0[i] = eqs0[i].subs(par0)
#     eqs0[i] = eqs0[i].subs(par1)
#     eqs0[i] = eqs0[i].subs(subc1)
#     eqs0[i] = eqs0[i].evalf(subs=par2)
#     i +=1
# # #eqs0 = [A[2,1],A[2,-1],B[2,1],B[2,-1]] # testline
# sol = solve(eqs0, unkn2)
# subc = sol
#
# subc = subc | subc1



subc = {}
for k in range(1,n+1):
    eqss = eqslist[(k-1)+n*(k-1)].copy() ##
    i = 0
    for ex in eqss:
        eqss[i] = ex.subs(par0)
        eqss[i] = eqss[i].subs(par1)
        eqss[i] = eqss[i].subs(subc)
        eqss[i] = eqss[i].evalf(subs=par2)
        i+=1
    sol0 = solve(eqss, unkn[k])
    subc = subc | sol0


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
#n = 1 # testline

Upsilon1p1 = Upsilon1p1.subs(z, x)
Upsilon1p1 = Upsilon1p1.subs(subc)

Upsilon1p1 = Upsilon1p1.subs(subf)

Upsilon1p1 = Upsilon1p1.subs(par0)
Upsilon1p1 = Upsilon1p1.subs(par1)
Upsilon1p1 = Upsilon1p1.subs(par2)
Upsilon1p1 = Upsilon1p1.evalf()  ##
#print(Upsilonp11)
#Upsilonp11 = Upsilonp11.subs(e**2, 0) # testline

phi1p1 = phi1p1.subs(z, x)
phi1p1 = phi1p1.subs(subc)

phi1p1 = phi1p1.subs(subf)

phi1p1 = phi1p1.subs(par0)
phi1p1 = phi1p1.subs(par1)
phi1p1 = phi1p1.subs(par2)
phi1p1 = phi1p1.evalf()

#phip11 = phip11.subs(e**2, 0) # testline

dphi1p1 = dphi1p1.subs(z, x)
dphi1p1 = dphi1p1.subs(subc)

dphi1p1 = dphi1p1.subs(subf)

dphi1p1 = dphi1p1.subs(par0)
dphi1p1 = dphi1p1.subs(par1)
dphi1p1 = dphi1p1.subs(par2)
dphi1p1 = dphi1p1.evalf()

#dphip11 = dphip11.subs(e**2, 0) # testline

expoa = 1 - 2 * I * e * df / (1 + I * e * df)  ##
expoa = expoa.series(e, n=n+1)
expoa = expoa.removeO()

h = (1 + e**2 * df**2)
hm = series(1/h, e ,n=n+1)  ##
hm = hm.removeO()

cr = e * d2f / h ** 3
cr = cr.series(e, n=n+1)
cr = cr.removeO()

w1 = symbols("w1", complex=True)

G1 = phi1p1+conjugate(phi1p1)-(Upsilon1p1+conjugate(phi1p1)-(w1-conjugate(w1))*conjugate(dphi1p1))*expoa


G1 = G1.subs(w1, x+I*e*f)


G1 = ser(G1, n+1)

G1 = G1.subs(subf)  ##
G1 = G1.subs(b, 2*pi)
#G1 = G1.subs(a, 10*10**(-9))

#G1 = G1.subs(sin(2*pi*x/a),0)
#G1 = G1.subs(cos(2*pi*x/a),0)
G1 = G1.subs(a, 10*10**(-9))

#G1 = G1.subs(b, 0)

#G1 = G1.subs(x, x*a)
#G1 = G1.subs(a, 10*10**(-9))

sigma_1nn = re(G1).evalf()  #eval
#print(sigma_1nn)

sigma_1tt_plus_sigma_1nn = re(4*phi1p1).evalf()

#print(sigma_1tt_plus_sigma_1nn)  ##  -0.I

sigma_1tt=sigma_1tt_plus_sigma_1nn-sigma_1nn

sigmatt = sigma_1tt

sigmatt = sigmatt.subs(e, 0.1)


scf = sigmatt.evalf(subs={x: 0})

plot(sigmatt, (x, -0.5, 0.5, 0.01))


x1 = np.arange(-0.5, 0.5, 0.01, dtype=float)
x2 = np.array([(10**(-9)*re(sigmatt)).subs(x, xi).evalf() for xi in x1])
#x2 = 10**(-9)*re(sigmatt).subs(x,x1)

# plt.plot(x1, x2)
# plt.show()

points = np.array([[x1[i], x2[i]] for i in range(len(x1)) ])

#plt.plot(points[:, 0], points[:, 1])
#plt.show()
out_file = open("points22","wb")
pickle.dump(points, out_file)
out_file.close()

#np.savetxt("points.csv", points, delimiter=",")