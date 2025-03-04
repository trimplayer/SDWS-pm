from sympy import *
import pickle
import numpy as np
import matplotlib.pyplot as plt



n_fur = 5

x, t, a, b, k, n, y_0 = symbols('x t a b k n y_0', real=True)
i = symbols('i', integer=True)
A = IndexedBase('A')
F = Symbol('F', real = True)

d_expr = im(cot(I * y_0)) + 1
w_expr = pi/a * x - I * y_0
f_expr = a/d_expr * (im(cot(w_expr)) - 1)

f_func = Lambda((x, a, y_0), f_expr)

fc_expr = (2/a) * Integral(f_expr.subs(x, t) * cos(2*k*pi*t/a), (t, -a/2, a/2))
fc_func = Lambda((a, y_0, k), fc_expr)

a_val = 1
y0_val = 0.5


coeffs = []
for j in range(1,n_fur+1):
    result_fc = fc_func(a_val, y0_val, j).evalf()
    coeffs.append(result_fc)
    #print(f"A_{j} = ", N(result_fc, 6))


F_expr = Sum(A[i] * a*cos(i*b*x), (i, 1, n_fur)).doit()
subs_dict = {A[i]: coeffs[i-1] for i in range(1, n_fur+1)}
#subs_dict.update({a: a_val})
result_F_expr = F_expr.subs(subs_dict)


dic_sub_pl = {b:2*pi,a:a_val}
F_x = Lambda((x), result_F_expr.subs(dic_sub_pl))

x_vals = np.linspace(-0.5, 0.5, 100)
f_values = [f_func(x, a_val, y0_val).evalf() for x in x_vals]
F_values = [F_x(x).evalf() for x in x_vals]

plt.plot(x_vals, f_values, 'r-', label=r'$f(x)$')
plt.plot(x_vals, F_values, 'b-', label=r'$F(x)$')
plt.xlabel('x', fontsize=14)
plt.ylabel('y', fontsize=14)
plt.legend(fontsize=14)
plt.grid(True)
plt.show()


print(result_F_expr,x)



out_file = open("fou_sub_f","wb")
pickle.dump(result_F_expr,out_file)
out_file.close()

out_file = open("fou_sub_df","wb")
pickle.dump(diff(result_F_expr,x),out_file)
out_file.close()

out_file = open("fou_sub_d2f","wb")
pickle.dump(diff(diff(result_F_expr,x),x),out_file)
out_file.close()


