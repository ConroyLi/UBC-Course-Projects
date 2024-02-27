from sympy import symbols, Eq, solve, latex, diff

# Define symbols
x1, x2, lam = symbols('x1 x2 lambda')

# Objective function
F = x1**2 + x2**2 + 10*x1 + 20*x2 + 25

# Constraint
constraint = x1 + x2

# Lagrangian
L = F + lam * (constraint - 0)  # Original constraint

# Partial derivatives
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

# Optimal values of x1, x2, lambda
opt_x1, opt_x2, opt_lam = solution[x1], solution[x2], solution[lam]
print('opt_x1=',latex(opt_x1))
print('opt_x2=',latex(opt_x2))
print('opt_lam=',latex(opt_lam))
# Calculate F_opt
F_opt = F.subs([(x1, opt_x1), (x2, opt_x2)])
print('F_opt=',latex(F_opt))

# For sensitivity analysis, change the constraint to x1 + x2 = 0.01
constraint_new = x1 + x2 - 0.01
L_new = F + lam * constraint_new

# Solve the new system of equations
solution_new = solve((diff(L_new, x1), diff(L_new, x2), diff(L_new, lam)), (x1, x2, lam))
print('L_new',latex(L_new))
print('partial_x1=',latex(diff(L_new, x1) ))
print('partial_x2=',latex(diff(L_new, x2)))
print('partial_lam=',latex(diff(L_new, lam)))
# New optimal values of x1, x2, lambda
opt_x1_new, opt_x2_new, opt_lam_new = solution_new[x1], solution_new[x2], solution_new[lam]
print('opt_x1=',latex(opt_x1_new))
print('opt_x2=',latex(opt_x2_new))
print('opt_lam=',latex(opt_lam_new))
# Calculate new F_opt
F_opt_new = F.subs([(x1, opt_x1_new), (x2, opt_x2_new)])
print('F_opt_new =',latex(F_opt_new ))
# Increase in F_opt
increase_in_F_opt = F_opt_new - F_opt
print('increase_in_F_opt =',latex(increase_in_F_opt ))
opt_x1, opt_x2, opt_lam, F_opt, opt_x1_new, opt_x2_new, opt_lam_new, F_opt_new, increase_in_F_opt
