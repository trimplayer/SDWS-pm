from sympy import *
in_file = open("rhs1.txt","r")
frexpr1 = in_file.read()
rhsn1 = parse_expr(frexpr1)
pprint(rhsn1)
print()