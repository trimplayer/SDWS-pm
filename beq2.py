from values import *
from sympy import *
from func import *
import pickle

n = 2

in_file = open("Upsilon1p1N","rb")
Upsilon1p1 = pickle.load(in_file)
in_file.close()

in_file = open("phi1p1N","rb")
phi1p1 = pickle.load(in_file)
in_file.close()

in_file = open("dphi1p1N","rb")
dphi1p1 = pickle.load(in_file)
in_file.close()



f = symbols('f', real=True)
df = symbols('df', real=True)
d2f = symbols('d2f', real=True)
e = symbols('e', real=True)
ka = symbols('kappa', real=True)
T = symbols('T', real=True)
mu = symbols('mu', real=True)
x = symbols('x', real=True)
y = symbols('y', real=True)
z = symbols('z')



subf = {f:sub_f, df:sub_df, d2f:sub_d2f}

#Upsilon1p1 = Upsilon1p1.subs({z: x-I*y, y:0})  # y = 0
qs0 = symbols('q_s0', real=True)
Upsilon1p1 = Upsilon1p1.subs(z, x)
Upsilon1p1 = Upsilon1p1.subs(qs0, 0)

phi1p1 = phi1p1.subs(z, x)
phi1p1 = phi1p1.subs(qs0, 0)

dphi1p1 = dphi1p1.subs(z, x)
dphi1p1 = dphi1p1.subs(qs0, 0)

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

expoa = 1 - 2 * I * e * df / (1 + I * e * df)  ##
expoa = expoa.series(e, n=n+1)
expoa = expoa.removeO()

h = (1 + e**2 * df**2)
hm = series((I + e * df)/h ,e ,n=n+1)
hm = hm.removeO()

cr = e * d2f / h ** 3
cr = cr.series(e, n=n+1)
cr = cr.removeO()

sigmass0 = symbols('sigma_ss0')
Ms = symbols("M_s")
sigmas0 = symbols("sigma_s0")

qs = (sigmass0*cr+Ms*cr*re(dus)+sigmas0*im(d2us*expoa))-I*(Ms*re(d2us*expoa)-sigmas0*cr*im(dus))
qs = ser(qs, n+1)


eqq4 = 2*mu*dus - (ka+1)*phi1p1 + qs
eqq4 = eqq4.expand()
eqq4 = collect(eqq4, e)


#qs0 = symbols('q_s0', real=True)  ##
#subsi = {dusij[0, 0]: ((ka+1)*T/4-qs0)/(2*mu), dusij[0,1]: 0, dusij[0,2]: 0}
subsi = {dusij[0, 0]: ((ka+1)*T/4)/(2*mu), dusij[0,1]: 0, dusij[0,2]: 0, dusij[0,3]: 0}

A = IndexedBase('A', real=True)
B = IndexedBase('B', real=True)
b = symbols('b', real=True)

sub_dus = (A[1,1] + I*A[1,(-1)])*cos(b*x) + (B[1,1] + I*B[1,(-1)])*sin(b*x)
sub_dus2 = (A[2,1] + I*A[2,(-1)])*cos(2*b*x) + (B[2,1] + I*B[2,(-1)])*sin(2*b*x)
subsi[dusij[1,0]] = sub_dus
subsi[dusij[1,1]] = diff(sub_dus, x)
subsi[dusij[1,2]] = diff(diff(sub_dus, x), x)
subsi[dusij[2,0]] = sub_dus2

#eqq41 = eqq4.coeff(e, 1)
eqq41 = eqq4.coeff(e, 1)
eqq42 = eqq4.coeff(e, 2)

eqq41 = eqq41.subs(subf)  ##
eqq41 = eqq41.subs(subsi)
#eqq41 = simplify(eqq41)
eqq41 = eqq41.expand()
eqq41 = eqq41.subs(exp(I*b*x), 1 / (exp(-I*b*x)))

eqq41 = eqq41.subs(exp(-I*b*x), cos(b*x)-I*sin(b*x))
eqq41 = eqq41.expand()
eqq41 = eqq41.collect([cos(b*x), sin(b*x)])

eqq42 = eqq42.subs(subf)  ##
eqq42 = eqq42.subs(subsi)
#eqq41 = simplify(eqq41)
eqq42 = eqq42.expand()
eqq42 = eqq42.subs(exp(I*b*x), 1 / (exp(-I*b*x)))
eqq42 = eqq42.subs(exp(2*I*b*x), 1 / (exp(-2*I*b*x)))

eqq42 = eqq42.subs(exp(-I*b*x), cos(b*x)-I*sin(b*x))
eqq42 = eqq42.subs(exp(-2*I*b*x), cos(2*b*x)-I*sin(2*b*x))
eqq42 = eqq42.expand()
eqq42 = eqq42.collect([cos(b*x), sin(b*x), cos(2*b*x), sin(2*b*x)])

cfc12 = eqq41.coeff(cos(b*x))
cfs12 = eqq41.coeff(sin(b*x))
eqs1 = [re(eqq41.coeff(cos(b*x))), im(eqq41.coeff(cos(b*x))), re(eqq41.coeff(sin(b*x))), im(eqq41.coeff(sin(b*x)))]
eqs12 = [re(eqq41.coeff(cos(2*b*x))), im(eqq41.coeff(cos(2*b*x))), re(eqq41.coeff(sin(2*b*x))), im(eqq41.coeff(sin(2*b*x)))]
eqs2 = [re(eqq42.coeff(cos(b*x))), im(eqq42.coeff(cos(b*x))), re(eqq42.coeff(sin(b*x))), im(eqq42.coeff(sin(b*x)))]
eqs22 = [re(eqq42.coeff(cos(2*b*x))), im(eqq42.coeff(cos(2*b*x))), re(eqq42.coeff(sin(2*b*x))), im(eqq42.coeff(sin(2*b*x)))]

eq1 = re(cfc12)
#eq2 = -1/2/mu*re(-(ka+1)*b*(-a*(sigmass0*mu+1/8*Ms*(T*ka+T-4*qs0))*b+(-a*T+(I*B[1, 1]+A[1, 1])*Ms+(I*A[1, -1]-B[1, -1])*sigmas0)*mu)+I*im(qs0)*a*b**2*sigmas0+re(qs0)*a*b**2*Ms-2*a*(sigmass0*mu+1/8*Ms*T*(ka+1))*b**2-4*(A[1, 1]+I*A[1, -1])*mu**2)
eq2 = -1/16/mu*re(-8*((-a*sigmass0*b+(I*B[1, 1]+A[1, 1])*Ms-T*a+sigmas0*(I*A[1, -1]-B[1, -1]))*mu-1/8*a*b*T*Ms*(ka+1))*(ka+1)*b-(I*32*A[1, -1]+32*A[1, 1])*mu**2-16*a*b**2*sigmass0*mu-2*a*b**2*Ms*T*(ka+1))
print(simplify(eq2-eq1))

eq3 = im(cfc12)
#eq4 = -1/2/mu*im(-(ka+1)*b*(-a*(sigmass0*mu+1/8*Ms*(T*ka+T-4*qs0))*b+(-a*T+(I*B[1, 1]+A[1, 1])*Ms+(I*A[1, -1]-B[1, -1])*sigmas0)*mu)+I*im(qs0)*a*b**2*sigmas0+re(qs0)*a*b**2*Ms-2*a*(sigmass0*mu+1/8*Ms*T*(ka+1))*b**2-4*(A[1, 1]+I*A[1, -1])*mu**2)
eq4 = -1/16/mu*im(-8*((-a*sigmass0*b+(I*B[1, 1]+A[1, 1])*Ms-T*a+sigmas0*(I*A[1, -1]-B[1, -1]))*mu-1/8*a*b*T*Ms*(ka+1))*(ka+1)*b-(I*32*A[1, -1]+32*A[1, 1])*mu**2-16*a*b**2*sigmass0*mu-2*a*b**2*Ms*T*(ka+1))
print(simplify(eq3-eq4))

eq5 = re(cfs12)
eq6 = 1/2/mu*re(-I*((-a*sigmass0*b+(I*B[1, 1]+A[1, 1])*Ms-T*a+sigmas0*(I*A[1, -1]-B[1, -1]))*mu-1/8*a*b*T*Ms*(ka+1))*(ka+1)*b+4*mu**2*(B[1, 1]+I*B[1, -1]))
print(simplify(eq5-eq6))

eq7 = im(cfs12)
eq8 = 1/2/mu*im(-I*((-a*sigmass0*b+(I*B[1, 1]+A[1, 1])*Ms-T*a+sigmas0*(I*A[1, -1]-B[1, -1]))*mu-1/8*a*b*T*Ms*(ka+1))*(ka+1)*b+4*mu**2*(B[1, 1]+I*B[1, -1]))
print(simplify(eq7-eq8))

out_file = open("cfc12","wb")
pickle.dump(cfc12, out_file)
out_file.close()

out_file = open("cfs12","wb")
pickle.dump(cfs12, out_file)
out_file.close()


eqslist = [eqs1, eqs12, eqs2, eqs22]
out_file = open("eqslist","wb")
pickle.dump(eqslist, out_file)
out_file.close()