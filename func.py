import sympy


def ser(expr, m):
    e = sympy.symbols('e')

    eq = sympy.expand(expr)
    eq = sympy.collect(eq,e)

    ret = sympy.Float(0)
    for j in range(m):
        ret += eq.coeff(e, j)*e**j

    return ret
