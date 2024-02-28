from sympy import symbols, diff, solve, Eq, latex

# Redefine symbols after reset
x, y, z, lam = symbols('x y z lambda')

# Objective function
F = (x - 2)**2 + y**2 + (z - 1)**2

# Constraint
constraint = z**2 - 4*x**2 - 2*y**2

# Lagrangian
L = F + lam * constraint

# Calculate partial derivatives
partial_x = diff(L, x)
partial_y = diff(L, y)
partial_z = diff(L, z)
partial_lam = diff(L, lam)

# Solve the system of equations
solution = solve((partial_x, partial_y, partial_z, partial_lam), (x, y, z, lam))

solution
print('F=',latex(F))
print('L=',latex(L))
print('partial_x1=',latex(partial_x ))
print('partial_x2=',latex(partial_y))
print('partial_x2=',latex(partial_z))
print('partial_lam=',latex(partial_lam))
print('solution=',latex(solution))