from values import *
from sympy import *
from func import *
import pickle

n = 1
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


