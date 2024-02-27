from sympy import symbols, pi, diff, solve, Eq, latex

# Define symbols
x1, x2, lam = symbols('x1 x2 lambda')

# Objective function (Negative of volume for minimization)
F = -pi * x1**2 * x2

# Constraint: Total surface area equals 24pi
A0 = 24*pi
constraint_eq = Eq(2*pi*x1**2 + 2*pi*x1*x2, A0)

# Lagrangian
L = F + lam * (2*pi*x1**2 + 2*pi*x1*x2 - A0)

# Calculate partial derivatives
partial_x1 = diff(L, x1)
partial_x2 = diff(L, x2)
partial_lam = diff(L, lam)
print('F=',latex(F))
print('L=',latex(L))
print('partial_x1=',latex(partial_x1 ))
print('partial_x2=',latex(partial_x2))
print('partial_lam=',latex(partial_lam))
# Solve the system of equations
solution = solve((partial_x1, partial_x2, partial_lam), (x1, x2, lam))

# Since there might be multiple solutions, filter out the physical ones (positive dimensions)
physical_solutions = [(x1_val, x2_val) for x1_val, x2_val, _ in solution if x1_val > 0 and x2_val > 0]

physical_solutions
print('physical_solutions=',latex(physical_solutions ))