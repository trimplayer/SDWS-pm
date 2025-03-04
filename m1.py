from values import *
from sympy import *
from func import *
import numpy as np
from matplotlib import pyplot as plt

#n = 3
nn = []
xx = []
xee = []
s_data = []

for n in range(1,3):
    ## 1
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

    ## 2

    phip1i1 = IndexedBase('Phip1_i_1')
    phip1ij1 = IndexedBase('Phi_p1_1')

    phip11 = Sum(e**i/factorial(i) * phip1i1[i], (i,0,n)).doit()
    for i in range(n+1):
        phip11 = phip11.subs(phip1i1[i], Sum((I*e*f)**j / factorial(j) * phip1ij1[i, j],(j,0,n)).doit())

    phip11 = ser(phip11, n+1)
    i, j = symbols('i j', cls=Idx)

    dphip1i1 = IndexedBase('dPhip1_i_1')

    dphip11 = Sum(e**i/factorial(i) * dphip1i1[i], (i,0,n)).doit()

    for i in range(n+1):
        dphip11 = dphip11.subs(dphip1i1[i], Sum((I*e*f)**j / factorial(j) * phip1ij1[i, j+1],(j,0,n)).doit())

    dphip11 = ser(dphip11, n+1)

    Upsilonp1i1 = IndexedBase('Upsilonp1_i_1')
    Upsilonp1ij1 = IndexedBase('Upsilon_p1_1')
    i, j = symbols('i j', cls=Idx)

    Upsilonp11 = Sum(e**i/factorial(i) * Upsilonp1i1[i], (i,0,n)).doit()

    for i in range(n+1):
        Upsilonp11 = Upsilonp11.subs(Upsilonp1i1[i], Sum((-I*e*f)**j / factorial(j) * Upsilonp1ij1[i, j],(j,0,n)).doit())

    Upsilonp11 = ser(Upsilonp11, n+1)

    expoa = 1 - 2*I*e*df/(1+I*e*df)
    expoa = expoa.series(e, n=n+1)
    expoa = expoa.removeO()

    w1 = symbols('w1')
    x = symbols('x', real=True)

    eqq1 = phip11 + conjugate(phip11) - (Upsilonp11 + conjugate(phip11) - (w1 - conjugate(w1)) * conjugate(dphip11)) * expoa - rhs1
    eqq1 = eqq1.subs(w1,x+I*e*f)

    eqq1 = ser(eqq1, n+1)



    ## 3

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
        eqq1l.append(eqq1.coeff(e, l))

    t = symbols('t', real=True)

    reqq1i = []

    for i in range(n):
        reqq1i.append(eqq1l[i+1].subs(subp))
        reqq1i[i] = reqq1i[i].subs(subf)
        reqq1i[i] = reqq1i[i].rewrite(exp(I*b*x))
        reqq1i[i] = reqq1i[i].expand()
        for j in range(1,n+1):
            reqq1i[i] = reqq1i[i].subs(exp(-j*I*b*x), 1/(t**j))
            reqq1i[i] = reqq1i[i].subs(exp(j*I*b*x), t**j)

        reqq1i[i] = collect(reqq1i[i], t)


    z = symbols('z')

    for i in range(1,n+1):
        Upsilonp11 = Upsilonp11.subs(Upsilonp1ij1[i, 0], reqq1i[i-1].coeff(t, i) * exp(i*I*b*z))
    Upsilonp11 = Upsilonp11.subs(subp)

    subphii = [0]
    for i in range(1,n+1):
        subphii.append(-1*(reqq1i[i-1].coeff(t, -i) * exp(-i*I*b*z)))

    for i in range(1,n+1):
        phip11 = phip11.subs(phip1ij1[i,0], subphii[i])
        subphidiff = subphii[i]
        for j in range(1,n+1):
            subphidiff = diff(subphidiff, z)
            phip11 = phip11.subs(phip1ij1[i, j], subphidiff)
    phip11 = phip11.subs(subp)

    for i in range(1,n+1):
        dphip11 = dphip11.subs(phip1ij1[i,0], subphii[i])
        subphidiff = subphii[i]
        for j in range(1,n): # n or n+1
            subphidiff = diff(subphidiff,z)
            dphip11 = dphip11.subs(phip1ij1[i, j], subphidiff)
    dphip11 = dphip11.subs(subp)




    ## 4

    qs0 = symbols('q_s0', real=True)
    Upsilonp11 = Upsilonp11.subs(z, x)
    Upsilonp11 = Upsilonp11.subs(qs0, 0)

    phip11 = phip11.subs(z, x)
    phip11 = phip11.subs(qs0, 0)

    dphip11 = dphip11.subs(z, x)
    dphip11 = dphip11.subs(qs0, 0)

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


    eqq4 = 2*mu*dus - (ka+1)*phip11 + qs
    eqq4 = eqq4.expand()
    eqq4 = collect(eqq4, e)

    A = IndexedBase('A', real=True)
    B = IndexedBase('B', real=True)
    b = symbols('b', real=True)

    subsi = {dusij[0, 0]: ((ka+1)*T/4)/(2*mu)}

    for k in range(1,n+2):
        subsi[dusij[0,k]] = 0

    for k in range(1,n+1):
        sub_dus = (A[k, 1] + I * A[k, (-1)]) * cos(k * b * x) + (B[k, 1] + I * B[k, (-1)]) * sin(k * b * x)
        subsi[dusij[k, 0]] = sub_dus
        sub_dus_d = diff(sub_dus,x)
        for l in range(1,n+1):

            subsi[dusij[k, l]] = sub_dus_d
            sub_dus_d = diff(sub_dus_d,x)

    eqq4l = [0]
    coslist = []
    for k in range(1,n+1):
        eqq4l.append(eqq4.coeff(e, k))
        eqq4l[k] = eqq4l[k].subs(subf)
        eqq4l[k] = eqq4l[k].subs(subsi)
        eqq4l[k] = eqq4l[k].expand()

        for l in range(n,0,-1):
            eqq4l[k] = eqq4l[k].subs(exp(l*I * b * x), 1 / (exp(-l*I * b * x)))
            eqq4l[k] = eqq4l[k].subs(exp(-l*I*b*x), cos(l*b*x)-I*sin(l*b*x))
        coslist.append(cos(k*b*x))
        coslist.append(sin(k*b*x))
        eqq4l[k] = eqq4l[k].expand()
        eqq4l[k] = eqq4l[k].collect(coslist)

    eqslist = []
    for k in range(1,n+1):
        for l in range(1,n+1):
            eqslist.append([re(eqq4l[k].coeff(cos(l*b*x))), im(eqq4l[k].coeff(cos(l*b*x))), re(eqq4l[k].coeff(sin(l*b*x))), im(eqq4l[k].coeff(sin(l*b*x)))])

    ## 5 stt
    unkn = [0]

    for k in range(1,n+1):
        unkn.append([A[k, 1],A[k, (-1)],B[k, 1],B[k, (-1)]])

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
    #par2 = {lam: 58.17 * 10 ** 9, mu: 26.13 * 10 ** 9, Ms: 6.0, a: 100 * 10 ** (-9), T: 0.1 * 10 ** 9, sigmass0: 0,
     #       sigmas0: 0}

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

    #print(subc)



    subf = {f:sub_f, df:sub_df, d2f:sub_d2f}

    z = symbols('z')


    Upsilonp11 = Upsilonp11.subs(z, x)
    Upsilonp11 = Upsilonp11.subs(subc)

    Upsilonp11 = Upsilonp11.subs(subf)

    Upsilonp11 = Upsilonp11.subs(par0)
    Upsilonp11 = Upsilonp11.subs(par1)
    Upsilonp11 = Upsilonp11.subs(par2)
    Upsilonp11 = Upsilonp11.evalf()


    phip11 = phip11.subs(z, x)
    phip11 = phip11.subs(subc)

    phip11 = phip11.subs(subf)

    phip11 = phip11.subs(par0)
    phip11 = phip11.subs(par1)
    phip11 = phip11.subs(par2)
    phip11 = phip11.evalf()



    dphip11 = dphip11.subs(z, x)
    dphip11 = dphip11.subs(subc)

    dphip11 = dphip11.subs(subf)

    dphip11 = dphip11.subs(par0)
    dphip11 = dphip11.subs(par1)
    dphip11 = dphip11.subs(par2)
    dphip11 = dphip11.evalf()

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

    G1 = phip11+conjugate(phip11)-(Upsilonp11+conjugate(phip11)-(w1-conjugate(w1))*conjugate(dphip11))*expoa


    G1 = G1.subs(w1, x+I*e*f)


    G1 = ser(G1, n+1)

    G1 = G1.subs(subf)  ##
    G1 = G1.subs(b, 2*pi)

    G1 = G1.subs(a, 10*10**(-9))



    sigma_1nn = re(G1).evalf()


    sigma_1tt_plus_sigma_1nn = re(4*phip11).evalf()



    sigma_1tt=sigma_1tt_plus_sigma_1nn-sigma_1nn

    #sigmatt = sigma_1tt
    sigmann = sigma_1nn
    xee = []
    xx=[]
    for ee in range(1,20):

    #sigmatt = sigmatt.subs(e, 0.1)
        #sigmatt = sigma_1tt
        #sigmatt = sigmatt.subs(e, ee*0.01)
        #scf = sigmatt.evalf(subs={x: 0})
        sigmann = sigma_1nn
        sigmann = sigmann.subs(e,ee*0.01)
        scf = sigmann.evalf(subs={x:0})
        ##scf = sigmatt.evalf(subs={x: 0})
        xee.append(ee*0.01)
        xx.append(scf)
    #scf = sigmatt.evalf(subs={x: 0})
    if s_data == []:
        s_data.append(xee)
    s_data.append(xx)
    #nn.append(n)
    #xx.append(scf)
    #plot(sigmatt, (x, -0.5, 0.5, 0.01))
    plt.plot(xee, xx,label='$n={n}$'.format(n=n))

#plt.plot(nn,xx)
#plt.plot(xee,xx)
s_data = np.array(s_data)
np.savetxt("snn_n15_a100_sig0.csv", s_data, delimiter=",")
#plt.plot(s_data[0, :], s_data[2, :])

plt.title('snn_n15')
plt.legend(loc='best')
plt.show()