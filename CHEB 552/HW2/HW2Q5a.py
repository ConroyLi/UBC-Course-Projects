from sympy import symbols, diff, sqrt, latex, Matrix, solve
# Redefine symbols for x1, x2, and introduce lambda1 and lambda2
x1, x2, lambda1, lambda2 = symbols('x1 x2 lambda1 lambda2')

# Define the objective function
F = x1**2 + x2**2 - 14*x1 - 6*x2 - 7

# Define the constraints
constraint1 = x1 + x2 - 2
constraint2 = x1 + 2*x2 - 3

# Form the Lagrangian
L = F + lambda1 * constraint1 + lambda2 * constraint2

# Compute the partial derivatives of the Lagrangian
partial_x1 = diff(L, x1)
partial_x2 = diff(L, x2)
partial_lambda1 = diff(L, lambda1)
partial_lambda2 = diff(L, lambda2)
second_partial_x1x1 = diff(partial_x1,x1)
second_partial_x1x2 = diff(partial_x1,x2)
second_partial_x2x2 = diff(partial_x2,x2)
# Solve the system of equations
solution = solve((partial_x1, partial_x2, partial_lambda1, partial_lambda2), (x1, x2, lambda1, lambda2))
H = Matrix([[second_partial_x1x1, second_partial_x1x2], 
            [second_partial_x1x2, second_partial_x2x2]])
F_min_value = F.subs([(x1, solution[x1]), (x2, solution[x2])])
# H_eva = H.subs([(x1, solution[x1]), (x2, solution[x2])])
H_eva = H.subs([(x1, 3), (x2, -1)])
eigenvalues = H_eva.eigenvals()
print('F=',latex(F))
print('L=',latex(L))
print('partial_x1=',latex(partial_x1 ))
print('partial_x2=',latex(partial_x2))
print('partial_lambda1=',latex(partial_lambda1))
print('partial_lambda2',latex(partial_lambda2))
print('solution=',latex(solution))
print('F_min_value=',latex(F_min_value))
print('H*=',latex(H_eva))
print('EV*=',latex(eigenvalues))