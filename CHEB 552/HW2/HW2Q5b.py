from sympy import symbols, diff, sqrt, latex, Matrix, solve
# Redefine symbols for lambda1 and lambda2 for the new problem
x1, x2, lambda1, lambda2 = symbols('x1 x2 lambda1 lambda2')
lambda1, lambda2 = symbols('lambda1 lambda2')

# Define the new objective function
F_new = x1**2 + 2*(x2 + 1)**2

# Define the new constraints
constraint1_new = -x1 + x2 - 2
constraint2_new = -x1 - x2 - 1

# Form the new Lagrangian
L_new = F_new + lambda1 * constraint1_new + lambda2 * constraint2_new

# Compute the partial derivatives of the new Lagrangian
partial_x1_new = diff(L_new, x1)
partial_x2_new = diff(L_new, x2)
partial_lambda1_new = diff(L_new, lambda1)
partial_lambda2_new = diff(L_new, lambda2)
second_partial_x1x1 = diff(partial_x1_new,x1)
second_partial_x1x2 = diff(partial_x1_new,x2)
second_partial_x2x2 = diff(partial_x2_new,x2)
# Solve the system of equations for the new problem
solution_new = solve((partial_x1_new, partial_x2_new, partial_lambda1_new, partial_lambda2_new), (x1, x2, lambda1, lambda2))
H = Matrix([[second_partial_x1x1, second_partial_x1x2], 
            [second_partial_x1x2, second_partial_x2x2]])
F_min_value = F_new.subs([(x1, solution_new[x1]), (x2, solution_new[x2])])
H_eva = H.subs([(x1, solution_new[x1]), (x2, solution_new[x2])])

eigenvalues = H_eva.eigenvals()
print('F=',latex(F_new))
print('L=',latex(L_new))
print('partial_x1=',latex(partial_x1_new ))
print('partial_x2=',latex(partial_x2_new))
print('partial_lambda1=',latex(partial_lambda1_new))
print('partial_lambda2',latex(partial_lambda2_new))
print('solution=',latex(solution_new))
print('F_min_value=',latex(F_min_value))
print('H*=',latex(H_eva))
print('EV*=',latex(eigenvalues))