from sympy import symbols, diff, sqrt, latex, Matrix, eye, trace,simplify

# Define symbolic variables
eta_xx, eta_xy, eta_yx, eta_yy, lambdas, mius = symbols('eta_xx eta_xy eta_yx eta_yy lambdas  mius')

F = Matrix([[1 + eta_xx, eta_yx],
            [eta_xy, 1 + eta_xy]])
FT = F.T
print(FT)
 
E = 0.5*(FT * F - eye(2))
print(E)

Pw = lambdas * E.trace() * F + 2 * mius * F * E
Pw = simplify(Pw)
print(latex(Pw))
