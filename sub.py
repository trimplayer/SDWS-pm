from values import *
from sympy import *
from func import ser

n = 1
x = symbols('x', real=True)
e = symbols('e', real=True)
f = symbols('f', real=True)
df = symbols('df', real=True)
d2f = symbols('d2f', real=True)
dict = {'e': e, 'f': f, 'df': df, 'd2f': d2f, 'x': x}

in_file = open("Upsilon1p1.txt.txt","r")
frexpr1 = in_file.read()
in_file.close()
Upsilon1p1 = parse_expr(frexpr1, local_dict=dict)

in_file = open("phi1p1.txt.txt","r")
frexpr1 = in_file.read()
in_file.close()
phi1p1 = parse_expr(frexpr1, local_dict=dict)

in_file = open("dphi1p1.txt.txt","r")
frexpr1 = in_file.read()
in_file.close()
dphi1p1 = parse_expr(frexpr1, local_dict=dict)

in_file = open("eqq11.txt.txt","r")
frexpr1 = in_file.read()
in_file.close()
eqq11 = parse_expr(frexpr1, local_dict=dict)

