from sympy import symbols, diff, sqrt, latex, Matrix, nsolve, log
# Redefine symbols for x1, x2, and introduce lambda1 and lambda2
x1, x2, lambda1, lambda2 = symbols('x1 x2 lambda1 lambda2', positive=True)

# Define the objective function
F = x1**2 + 1.5*x2**2 - 4*x1 - 7*x2 +x1*x2 + 9 - log(x1) - log(x2)

# Define the constraints
constraint1 = 4 - x1*x2
constraint2 = 2*x1 - x2 

# Form the Lagrangian
La = F 
Lb = F + lambda1 * constraint1 
Lc = F + lambda1 * constraint1 + lambda2 * constraint2
# Compute the partial derivatives of the Lagrangian
partial_x1_a = diff(La, x1)
partial_x2_a = diff(La, x2)

partial_x1_b = diff(Lb, x1)
partial_x2_b = diff(Lb, x2)
partial_lambda1_b = diff(Lb, lambda1)

partial_x1_c = diff(Lc, x1)
partial_x2_c = diff(Lc, x2)
partial_lambda1_c = diff(Lc, lambda1)
partial_lambda2_c = diff(Lc, lambda2)

second_partial_x1x1 = diff(partial_x1_a,x1)
second_partial_x1x2 = diff(partial_x1_a,x2)
second_partial_x2x2 = diff(partial_x2_a,x2)
# nsolve the system of equations
solution_a = nsolve((partial_x1_a, partial_x2_a,), (x1, x2), (1,1))
solution_b = nsolve((partial_x1_b, partial_x2_b, partial_lambda1_b), (x1, x2, lambda1),(1,1,1))
solution_c = nsolve((partial_x1_c, partial_x2_c, partial_lambda1_c, partial_lambda2_c), (x1, x2, lambda1, lambda2),(1,1,1,1))

H = Matrix([[second_partial_x1x1, second_partial_x1x2], 
            [second_partial_x1x2, second_partial_x2x2]])
F_min_value_a = F.subs([(x1, solution_a[0]), (x2, solution_a[1])])
F_min_value_b = F.subs([(x1, solution_b[0]), (x2, solution_b[1])])
F_min_value_c = F.subs([(x1, solution_c[0]), (x2, solution_c[1])])

H_a_eva = H.subs([(x1, solution_a[0]), (x2, solution_a[1])])
H_b_eva = H.subs([(x1, solution_b[0]), (x2, solution_b[1])])
H_c_eva = H.subs([(x1, solution_c[0]), (x2, solution_c[1])])
eigenvalues_a = H_a_eva.eigenvals()
eigenvalues_b = H_b_eva.eigenvals()
eigenvalues_c = H_c_eva.eigenvals()
print('F=',latex(F))
print('La=',latex(La))
print('Lb=',latex(Lb))
print('Lc=',latex(Lc))

print('partial_x1_a=',latex(partial_x1_a))
print('partial_x2_a',latex(partial_x2_a))

print('partial_x1_b=',latex(partial_x1_b))
print('partial_x2_b',latex(partial_x2_b))
print('partial_lambda1_b=',latex(partial_lambda1_b))

print('partial_x1_c=',latex(partial_x1_c))
print('partial_x2_c',latex(partial_x2_c))
print('partial_lambda1_c=',latex(partial_lambda1_c))
print('partial_lambda2_c=',latex(partial_lambda2_c))
print('solutiona=',latex(solution_a))
print('solutionb=',latex(solution_b))
print('solutionc=',latex(solution_c))
print('F_min_value_a=',latex(F_min_value_a))
print('F_min_value_b=',latex(F_min_value_b))
print('F_min_value_c=',latex(F_min_value_c))
print('Ha*=',latex(H_a_eva))
print('EVa*=',latex(eigenvalues_a))
print('Hb*=',latex(H_b_eva))
print('EVb*=',latex(eigenvalues_b))
print('Hc*=',latex(H_c_eva))
print('EVc*=',latex(eigenvalues_c))