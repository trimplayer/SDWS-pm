from values import *
from sympy import *
from func import *
import pickle

n = 4
x = symbols('x', real=True)
e = symbols('e', real=True)
f = symbols('f', real=True)
df = symbols('df', real=True)
d2f = symbols('d2f', real=True)


in_file = open("Upsilon1p1","rb")
Upsilon1p1 = pickle.load(in_file)
in_file.close()


in_file = open("phi1p1","rb")
phi1p1 = pickle.load(in_file)
in_file.close()


in_file = open("dphi1p1","rb")
dphi1p1 = pickle.load(in_file)
in_file.close()


in_file = open("eqq11","rb")
eqq11 = pickle.load(in_file)
in_file.close()

subf = {f:sub_f, df:sub_df, d2f:sub_d2f}

ka = symbols('kappa', real=True)
T = symbols('T', real=True)
mu = symbols('mu', real=True)
qs0 = symbols('q_s0', real=True)

dusij = IndexedBase("du_s")
Upsilonp1ij1 = IndexedBase('Upsilon_p1_1')
phip1ij1 = IndexedBase('Phi_p1_1')

subp = {}
subp[dusij[0, 0]] = ((ka+1)*T/4 - qs0) / (2*mu)
subp[phip1ij1[0, 0]] = T/4
subp[Upsilonp1ij1[0, 0]] = T/4

for l in range(1,n+2):
    subp[dusij[0, l]] = 0
    subp[phip1ij1[0, l]] = 0
    subp[Upsilonp1ij1[0, l]] = 0


for k in range(1,n+1):
    A = IndexedBase('A', real=True)
    B = IndexedBase('B', real=True)

    b = symbols('b', real=True)

    sub_dus = (A[k,1] + I*A[k,(-1)])*cos(k*b*x) + (B[k,1] + I*B[k,(-1)])*sin(k*b*x)

    subp[dusij[k,0]] = sub_dus

    subp_dusd = sub_dus
    for l in range(1,n+1):
        subp_dusd = subp_dusd.diff(x)
        subp[dusij[k,l]] = subp_dusd
    subp[phip1ij1[k,0]] = 0  ##
    subp[Upsilonp1ij1[k, 0]] = 0
    for l in range(1, n + 1):
        subp[phip1ij1[k, l]] = 0
        subp[Upsilonp1ij1[k, l]] = 0





eqq1l =[]
for l in range(n+1):
    eqq1l.append(eqq11.coeff(e, l))

# reqq11 = eqq11.subs(subp)
# reqq11 = reqq11.subs(subf)
# #reqq11 = simplify(reqq11)
# reqq11 = reqq11.rewrite(exp(I*b*x))

#reqq11 = eqq1l[1].subs(subp) # todo: for any power
#reqq11 = reqq11.subs(subf)
#reqq11 = reqq11.rewrite(exp(I*b*x))
#reqq12 = eqq1l[2].subs(subp)
#reqq12 = reqq12.subs(subf)
#reqq12 = reqq12.rewrite(exp(I*b*x))

t = symbols('t', real=True)

#reqq11 = reqq11.expand()
#reqq11 = reqq11.subs(exp(-I*b*x), 1/t)
#reqq11 = reqq11.subs(exp(I*b*x), t)
#reqq11 = reqq11.subs(exp(-2*I*b*x), 1/(t**2)) # todo: for any power
#reqq11 = reqq11.subs(exp(2*I*b*x), t**2)
#reqq11 = collect(reqq11, t)
#reqq12 = reqq12.expand()
#reqq12 = reqq12.subs(exp(-I*b*x), 1/t)
#reqq12 = reqq12.subs(exp(I*b*x), t)
#reqq12 = reqq12.subs(exp(-2*I*b*x), 1/(t**2))
#reqq12 = reqq12.subs(exp(2*I*b*x), t**2)
#reqq12 = collect(reqq12, t)

reqq1i = []

for i in range(n):
    reqq1i.append(eqq1l[i+1].subs(subp))
    reqq1i[i] = reqq1i[i].subs(subf)
    reqq1i[i] = reqq1i[i].rewrite(exp(I*b*x))
    reqq1i[i] = reqq1i[i].expand()
    for j in range(1,n+1):
        reqq1i[i] = reqq1i[i].subs(exp(-j*I*b*x), 1/(t**j))
        reqq1i[i] = reqq1i[i].subs(exp(j*I*b*x), t**j)

    reqq1i[i] = collect(reqq1i[i], t) ##coolect?

#  upsilon [1, 0]
z = symbols('z')
#Upsilonp11 = Upsilonp11.subs(Upsilonp1ij1[1, 0], reqq1i[0].coeff(t, 1) * exp(I*b*z)) # todo: for any power
#Upsilonp11 = Upsilonp11.subs(Upsilonp1ij1[2, 0], reqq1i[1].coeff(t, 2) * exp(2*I*b*z))
#Upsilonp11 = Upsilonp11.subs(subp)

for i in range(1,n+1):
    Upsilon1p1 = Upsilon1p1.subs(Upsilonp1ij1[i, 0], reqq1i[i-1].coeff(t, i) * exp(i*I*b*z))
Upsilon1p1 = Upsilon1p1.subs(subp)

out_file = open("Upsilon1p1N","wb")
pickle.dump(Upsilon1p1, out_file)
out_file.close()


#subphi = -1*(reqq11.coeff(t, -1) * exp(-I*b*z)) # todo: for any power
#subphi2 = -1*(reqq12.coeff(t, -2) * exp(-2*I*b*z))

subphii = [0]
for i in range(1,n+1):
    subphii.append(-1*(reqq1i[i-1].coeff(t, -i) * exp(-i*I*b*z)))

#phip11 = phip11.subs(phip1ij1[1,0], subphi)
#phip11 = phip11.subs(phip1ij1[2,0], subphi2)
#phip11 = phip11.subs(phip1ij1[1,1], diff(subphi, z))
#phip11 = phip11.subs(phip1ij1[1,2], diff(subphi2, z))
#phip11 = phip11.subs(subp)

for i in range(1,n+1):
    phi1p1 = phi1p1.subs(phip1ij1[i,0], subphii[i])
    subphidiff = subphii[i]
    for j in range(1,n+1):
        subphidiff = diff(subphidiff, z)
        phi1p1 = phi1p1.subs(phip1ij1[i, j], subphidiff)
phi1p1 = phi1p1.subs(subp)

out_file = open("phi1p1N","wb")
pickle.dump(phi1p1,out_file)
out_file.close()

#dphip11 = dphip11.subs(phip1ij1[1,0], subphi)  # todo: for any power
#dphip11 = dphip11.subs(phip1ij1[1,1], diff(subphi, z))
#dphip11 = dphip11.subs(phip1ij1[2,0], subphi2)
#dphip11 = dphip11.subs(phip1ij1[2,1], diff(subphi2, z))
#dphip11 = dphip11.subs(subp)

for i in range(1,n+1):
    dphi1p1 = dphi1p1.subs(phip1ij1[i,0], subphii[i])
    subphidiff = subphii[i]
    for j in range(1,n): # n or n+1
        subphidiff = diff(subphidiff,z)
        dphi1p1 = dphi1p1.subs(phip1ij1[i, j], subphidiff)
dphi1p1 = dphi1p1.subs(subp)

out_file = open("dphi1p1N","wb")
pickle.dump(dphi1p1, out_file)
out_file.close()