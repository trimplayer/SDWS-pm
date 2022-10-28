from sympy import *
x, a, b = symbols('x a b')
sub_f = -a * cos(b * x)
sub_df = a*b*sin(b*x)
sub_d2f = a*b**2*cos(b*x)
