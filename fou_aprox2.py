from values import *
from sympy import *
from func import *
import numpy as np
import pickle
#from fourier_fun import n_fur
from matplotlib import pyplot as plt

n_fur = 5 # fourier order              ######## don't forget to change
n_max = 1  # max approximation order
s_data = [] # data to save in file


# fourier
in_file = open("fou_sub_f","rb")
fou_sub_f = pickle.load(in_file)
in_file.close()
in_file = open("fou_sub_df","rb")
fou_sub_df = pickle.load(in_file)
in_file.close()
in_file = open("fou_sub_d2f","rb")
fou_sub_d2f = pickle.load(in_file)
in_file.close()

for n in range(1, n_max+1):

    #
    # 1 - rhs
    #

    # symbols and indexes
    e = symbols('e', real=True)
    f = symbols('f', real=True)
    df = symbols('df', real=True)
    d2f = symbols('d2f', real=True)
    sigmass0 = symbols('sigma_ss0')
    sigmas0 = symbols("sigma_s0")
    Ms = symbols("M_s")
    w1 = symbols('w1')  # temporary symbol for substitution
    x = symbols('x', real=True)
    ka = symbols('kappa', real=True)
    T = symbols('T', real=True)
    mu = symbols('mu', real=True)
    qs0 = symbols('q_s0', real=True)
    t = symbols('t', real=True)
    z = symbols('z')
    nu = symbols("nu")
    lam = symbols("lambda")

    dusi = IndexedBase("dus_i")
    dusij = IndexedBase("du_s")
    d2usi = IndexedBase("dus_i")
    phip1i1 = IndexedBase('Phip1_i_1')
    phip1ij1 = IndexedBase('Phi_p1_1')
    dphip1i1 = IndexedBase('dPhip1_i_1')
    Upsilonp1i1 = IndexedBase('Upsilonp1_i_1')
    Upsilonp1ij1 = IndexedBase('Upsilon_p1_1')
    A = IndexedBase('A', real=True)
    B = IndexedBase('B', real=True)
    b = symbols('b', real=True)
    i, j = symbols('i j', cls=Idx)

    #subf = {f: sub_f, df: sub_df, d2f: sub_d2f}
    subf = {f: fou_sub_f, df: fou_sub_df, d2f: fou_sub_d2f}
    # series

    dus = Sum((e ** i / factorial(i)) * dusi[i], (i, 0, n)).doit()

    for k in range(n + 1):
        dus = dus.subs(dusi[k], Sum(((I * e * f) ** j / factorial(j)) * dusij[k, j], (j, 0, n)).doit())

    dus = ser(dus, n + 1)

    d2us = Sum((e ** i / factorial(i)) * d2usi[i], (i, 0, n)).doit()

    for k in range(n + 1):
        d2us = d2us.subs(d2usi[k], Sum(((I * e * f) ** j / factorial(j)) * dusij[k, j + 1], (j, 0, n)).doit())

    d2us = ser(d2us, n + 1)

    expoa = (1 + I * e * df) / sqrt(1 + e ** 2 * df ** 2)
    expoa = expoa.series(e, n=n + 1)
    expoa = expoa.removeO()

    h = (1 + e ** 2 * df ** 2)  #
    hm = series((I + e * df) / h, e, n=n + 1)
    hm = hm.removeO()

    cr = e * d2f / h ** 3
    cr = cr.series(e, n=n + 1)
    cr = cr.removeO()

    # rhs1 - qs function of model
    rhs1 = (sigmass0 * cr + Ms * cr * re(dus) + sigmas0 * im(d2us * expoa)) - I * (
            Ms * re(d2us * expoa) - sigmas0 * cr * im(dus))
    rhs1 = ser(rhs1, n + 1)

    #
    # 2 - eq with phi and upsilon functions
    #

    # small parameter epsilon and Taylor series for phi function

    i, j = symbols('i j', cls=Idx)
    phip11 = Sum(e ** i / factorial(i) * phip1i1[i], (i, 0, n)).doit()

    for i in range(n + 1):
        phip11 = phip11.subs(phip1i1[i], Sum((I * e * f) ** j / factorial(j) * phip1ij1[i, j], (j, 0, n)).doit())

    phip11 = ser(phip11, n + 1)

    # dPhi series

    i, j = symbols('i j', cls=Idx)
    dphip11 = Sum(e ** i / factorial(i) * dphip1i1[i], (i, 0, n)).doit()

    for i in range(n + 1):
        dphip11 = dphip11.subs(dphip1i1[i], Sum((I * e * f) ** j / factorial(j) * phip1ij1[i, j + 1], (j, 0, n)).doit())

    dphip11 = ser(dphip11, n + 1)

    # upsilon function series

    i, j = symbols('i j', cls=Idx)
    Upsilonp11 = Sum(e ** i / factorial(i) * Upsilonp1i1[i], (i, 0, n)).doit()

    for i in range(n + 1):
        Upsilonp11 = Upsilonp11.subs(Upsilonp1i1[i],
                                     Sum((-I * e * f) ** j / factorial(j) * Upsilonp1ij1[i, j], (j, 0, n)).doit())

    Upsilonp11 = ser(Upsilonp11, n + 1)

    # exp(alpha) with minus

    expoa = 1 - 2 * I * e * df / (1 + I * e * df)
    expoa = expoa.series(e, n=n + 1)
    expoa = expoa.removeO()

    # equation of qs

    eqq1 = phip11 + conjugate(phip11) - (
            Upsilonp11 + conjugate(phip11) - (w1 - conjugate(w1)) * conjugate(dphip11)) * expoa - rhs1
    eqq1 = eqq1.subs(w1, x + I * e * f)

    eqq1 = ser(eqq1, n + 1)                         ####### checked with hands

    #
    # 3 - substitution
    #

    subp = {dusij[0, 0]: ((ka + 1) * T / 4 - qs0) / (2 * mu), phip1ij1[0, 0]: T / 4, Upsilonp1ij1[0, 0]: T / 4}

    for l in range(1, n + 2):
        subp[dusij[0, l]] = 0
        subp[phip1ij1[0, l]] = 0
        subp[Upsilonp1ij1[0, l]] = 0

    # forming substitution

    for k in range(1, n + 1):
        #sub_dus = (A[k, 1] + I * A[k, (-1)]) * cos(k * b * x) + (B[k, 1] + I * B[k, (-1)]) * sin(k * b * x) # for non-fourier series
        sub_dus = 0
        for k_f in range (1,n_fur+1):
            sub_dus += (A[k, k_f] + I * A[k, (-k_f)]) * cos(k_f * b * x) + (B[k, k_f] + I * B[k, (-k_f)]) * sin(k_f * b * x)


        subp[dusij[k, 0]] = sub_dus
        subp_dusd = sub_dus
        for l in range(1, n + 1):
            subp_dusd = subp_dusd.diff(x)
            subp[dusij[k, l]] = subp_dusd
        subp[phip1ij1[k, 0]] = 0
        subp[Upsilonp1ij1[k, 0]] = 0
        for l in range(1, n + 1):
            subp[phip1ij1[k, l]] = 0
            subp[Upsilonp1ij1[k, l]] = 0

    #

    eqq1l = []
    for l in range(n + 1):
        eqq1l.append(eqq1.coeff(e, l))

    reqq1i = []



    for i in range(n):
        reqq1i.append(eqq1l[i+1].subs(subp))
        reqq1i[i] = reqq1i[i].subs(subf)
        reqq1i[i] = reqq1i[i].rewrite(exp(I * b * x))
        reqq1i[i] = reqq1i[i].expand()

        for j in range(n_fur,0,-1):
            reqq1i[i] = reqq1i[i].subs(exp(-j * I * b * x), 1 / (t ** j))
            reqq1i[i] = reqq1i[i].subs(exp(j * I * b * x), t ** j)

        reqq1i[i] = collect(reqq1i[i], t)

    # potentials substitution

    subupsilon = [0]
    for k in range(1,n+1):
        subupsilon.append(0)
        for k_f in range(1,n_fur+1):

            subupsilon[k] += (reqq1i[k-1].coeff(t,k_f)*exp(k_f * I * b * z))



    for i in range(1, n + 1):
        Upsilonp11 = Upsilonp11.subs(Upsilonp1ij1[i, 0], subupsilon[i])
        subupsilondiff = subupsilon[i]
        for j in range(1,n + 1):
            subupsilondiff = diff(subupsilondiff, z)
            Upsilonp11 = Upsilonp11.subs(Upsilonp1ij1[i, j], subupsilondiff)
    Upsilonp11 = Upsilonp11.subs(subp)              #####   checked manually

    subphii = [0]
    for k in range(1,n+1):
        subphii.append(0)
        for k_f in range(1,n_fur+1):

            subphii[k] += (-1*(reqq1i[k-1].coeff(t,-k_f)*exp(-k_f*I*b*z)))

    for i in range(1, n + 1):
        phip11 = phip11.subs(phip1ij1[i, 0], subphii[i])
        subphidiff = subphii[i]
        for j in range(1, n + 1):
            subphidiff = diff(subphidiff, z)
            phip11 = phip11.subs(phip1ij1[i, j], subphidiff)
    phip11 = phip11.subs(subp)                                      #######  checked manually

    for i in range(1, n + 1):
        dphip11 = dphip11.subs(phip1ij1[i, 0], subphii[i])
        subphidiff = subphii[i]
        for j in range(1, n + 1):  # n or n+1
            subphidiff = diff(subphidiff, z)
            dphip11 = dphip11.subs(phip1ij1[i, j], subphidiff)
    dphip11 = dphip11.subs(subp)

    #
    # 4 - rewrite equasions
    #

    Upsilonp11 = Upsilonp11.subs(z, x)  # substitutions and series
    Upsilonp11 = Upsilonp11.subs(qs0, 0)

    phip11 = phip11.subs(z, x)
    phip11 = phip11.subs(qs0, 0)

    dphip11 = dphip11.subs(z, x)
    dphip11 = dphip11.subs(qs0, 0)

    i, j = symbols('i j', cls=Idx)
    dus = Sum((e ** i / factorial(i)) * dusi[i], (i, 0, n)).doit()
    for k in range(n + 1):
        dus = dus.subs(dusi[k], Sum(((I * e * f) ** j / factorial(j)) * dusij[k, j], (j, 0, n)).doit())
    dus = ser(dus, n + 1)

    d2us = Sum((e ** i / factorial(i)) * d2usi[i], (i, 0, n)).doit()
    for k in range(n + 1):
        d2us = d2us.subs(d2usi[k], Sum(((I * e * f) ** j / factorial(j)) * dusij[k, j + 1], (j, 0, n)).doit())
    d2us = ser(d2us, n + 1)

    expoa = 1 - 2 * I * e * df / (1 + I * e * df)
    expoa = expoa.series(e, n=n + 1)
    expoa = expoa.removeO()
    h = (1 + e ** 2 * df ** 2)
    hm = series((I + e * df) / h, e, n=n + 1)
    hm = hm.removeO()
    cr = e * d2f / h ** 3
    cr = cr.series(e, n=n + 1)
    cr = cr.removeO()

    # eq
    qs = (sigmass0 * cr + Ms * cr * re(dus) + sigmas0 * im(d2us * expoa)) - I * (
                Ms * re(d2us * expoa) - sigmas0 * cr * im(dus))

    qs = ser(qs, n + 1)

    eqq4 = 2 * mu * dus - (ka + 1) * phip11 + qs
    eqq4 = eqq4.expand()
    eqq4 = collect(eqq4, e)         ###### checked manually


    # forming subs
    subsi = {dusij[0, 0]: ((ka + 1) * T / 4) / (2 * mu)}

    for k in range(1, n + 2):
        subsi[dusij[0, k]] = 0

    for k in range(1, n + 1):
        #sub_dus = (A[k, 1] + I * A[k, (-1)]) * cos(k * b * x) + (B[k, 1] + I * B[k, (-1)]) * sin(k * b * x) # origin for non-fourier
        sub_dus = 0
        for k_f in range(1,n_fur+1):
            sub_dus += (A[k, k_f] + I * A[k, (-k_f)]) * cos(k_f * b * x) + (B[k, k_f] + I * B[k, (-k_f)]) * sin(k_f * b * x)
        subsi[dusij[k, 0]] = sub_dus
        sub_dus_d = diff(sub_dus, x)
        for l in range(1, n + 1):
            subsi[dusij[k, l]] = sub_dus_d
            sub_dus_d = diff(sub_dus_d, x)

    # get equations to solve on next step
    eqq4l = [0]
    #coslist = []


    for k in range(1, n + 1):

        eqq4l.append(eqq4.coeff(e, k))
        eqq4l[k] = eqq4l[k].subs(subf)
        eqq4l[k] = eqq4l[k].subs(subsi)
        eqq4l[k] = eqq4l[k].expand()
        eqq4l[k] = eqq4l[k].rewrite(exp(I * b * x))

        trig_sub = {}
        for k_f in range(n_fur,0,-1):
           trig_sub[exp(-k_f * I * b * x)] = 1 / (t ** k_f)
           trig_sub[exp(k_f * I * b * x)] = t ** k_f
        eqq4l[k] = eqq4l[k].subs(trig_sub)
        eqq4l[k] = eqq4l[k].expand()
        eqq4l[k] = eqq4l[k].evalf()
        eqq4l[k] = eqq4l[k].collect(t)







    eqslist = []

    for k in range(1,n+1):
        cfc = [0]
        cfs = [0]
        for k_f in range(1,n_fur+1):
            cfc.append(eqq4l[k].coeff(t,k_f))
            cfs.append(eqq4l[k].coeff(t, -k_f))


        for k_f in range(1,n_fur+1):
            eqslist.append([re(cfc[k_f]),im(cfc[k_f]),re(cfs[k_f]),im(cfs[k_f])])

    #
    # 5 - get coeffs and stress components
    #

    unkn = [] # forming unknown coeffs matrix

    for k in range(1, n + 1):
        #unkn.append([A[k, 1], A[k, (-1)], B[k, 1], B[k, (-1)]]) # origin for non-fourier

        for k_f in range(1,n_fur+1):

            unkn.append(A[k, k_f])
            unkn.append(A[k, (-k_f)])
            unkn.append(B[k, k_f])
            unkn.append(B[k, (-k_f)])

    par0 = {qs0: 0}

    # first subs parameters
    par1 = {b: 2 * pi / a, ka: 3 - 4 * (lam / (2 * (lam + mu))), x: x * a, nu: lam / (2 * (lam + mu))}

    # second subs - Ms; a; T; sigma
    par2 = {lam: 58.17 * 10 ** 9, mu: 26.13 * 10 ** 9, Ms: 6.099, a: 10 * 10 ** (-9), T: 0.1 * 10 ** 9, sigmass0: 1,
            sigmas0: 1}  # origin

    # solving equations

    subc = {}  # coeffs

    for kn in range(1,n+1):
        eqss = []
        for k in range(1,n_fur+1):
            for kir in range(4):
                cure = eqslist[(k - 1)+(kn-1)*n_fur][kir]
                cure = cure.subs(par0)
                cure = cure.subs(par1)
                cure = cure.subs(subc)
                cure = cure.evalf(subs=par2)
                eqss.append(cure)

        sol0 = solve(eqss, unkn)

        subc = subc | sol0

    #subf = {f: sub_f, df: sub_df, d2f: sub_d2f}  # function and derivatives
    subf = {f: fou_sub_f, df: fou_sub_df, d2f: fou_sub_d2f}

    # subs into potentials
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


    expoa = 1 - 2 * I * e * df / (1 + I * e * df)  # get series
    expoa = expoa.series(e, n=n + 1)
    expoa = expoa.removeO()
    h = (1 + e ** 2 * df ** 2)
    hm = series(1 / h, e, n=n + 1)
    hm = hm.removeO()
    cr = e * d2f / h ** 3
    cr = cr.series(e, n=n + 1)
    cr = cr.removeO()

    # stress
    G1 = phip11 + conjugate(phip11) - (
                Upsilonp11 + conjugate(phip11) - (w1 - conjugate(w1)) * conjugate(dphip11)) * expoa

    G1 = G1.subs(w1, x + I * e * f)
    G1 = ser(G1, n + 1)
    #G1 = G1.subs(subf)
    #G1 = G1.subs(b, par1[b])
    #G1 = G1.subs(a, par2[a])


    sigma_1nn = re(G1).evalf()
    sigma_1tt_plus_sigma_1nn = re(4 * phip11).evalf()
    sigma_1tt = sigma_1tt_plus_sigma_1nn - sigma_1nn

    # G1 -> stress

    sigmatt = sigma_1tt
    sigmann = sigma_1nn

    scf = sigmatt.evalf(subs={x: 0})

    sigmatt = sigmatt.subs(e, 0.1)
    #sigmatt = sigmatt.subs(f, 0)
    #sigmatt = sigmatt.subs(df,0)

    sigmatt = sigmatt.subs(f, fou_sub_f)
    sigmatt = sigmatt.subs(df,fou_sub_df)
    sigmatt = sigmatt.subs(b, par1[b])
    sigmatt = sigmatt.subs(a, par2[a])
    #sigmatt = sigmatt.subs(a, -1)







    plot(sigmatt, (x, -0.5, 0.5))

    # data saving and plotting
    xee = []
    xx = []
    for ee in range(1, 11):
        # block for sigma tt


        sigmatt = sigma_1tt
        sigmatt = sigmatt.subs(f, fou_sub_f)
        sigmatt = sigmatt.subs(df, fou_sub_df)
        sigmatt = sigmatt.subs(b, par1[b])
        sigmatt = sigmatt.subs(a, par2[a])
        sigmatt = sigmatt.subs(e, ee*0.01)
        scf = sigmatt.evalf(subs={x: 0})

        # block for sigma nn

        #sigmann = sigma_1nn
        #sigmann = sigmann.subs(e, ee * 0.01)
        #scf = sigmann.evalf(subs={x: 0})


        xee.append(ee * 0.01)
        xx.append(scf)

    if s_data == []:
        s_data.append(xee)
    s_data.append(xx)
    #plt.plot(xee, xx, label='$n={n}$'.format(n=n))

# data saving
s_data = np.array(s_data)
#np.savetxt("stt_fou_y05_n1.csv", s_data, delimiter=",")


#plt.title('snn_n15')
#plt.legend(loc='best')
plt.show()


