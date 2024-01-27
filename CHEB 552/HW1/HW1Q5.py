from sympy import symbols, diff, solve, Matrix, latex, N

# Define symbols
x1, x2, x3 = symbols('x1 x2 x3')

# Function a
f_a = x1**4 + 12*x2**3 - 15*x1**2 - 56*x2 + 60
print('F_a=',latex(f_a))
# Partial derivatives for Function a
partial_f_a_x1 = diff(f_a, x1)
partial_f_a_x2 = diff(f_a, x2)
print('dF_a1=',latex(partial_f_a_x1))
print('dF_a2=',latex(partial_f_a_x2))
# Stationary points for Function a
stationary_points_a = solve((partial_f_a_x1, partial_f_a_x2), (x1, x2))
print('sta_a=',latex(stationary_points_a))

# Hessian matrix for Function a
second_partial_f_a_x1x1 = diff(partial_f_a_x1, x1)
second_partial_f_a_x1x2 = diff(partial_f_a_x1, x2)
second_partial_f_a_x2x2 = diff(partial_f_a_x2, x2)
Hessian_a = Matrix([
    [second_partial_f_a_x1x1, second_partial_f_a_x1x2],
    [second_partial_f_a_x1x2, second_partial_f_a_x2x2]
])
print('H_a=',latex(Hessian_a))
# Function b
f_b = x1**2 + x2**2 + x3**2 - 4*x1*x2
print('F_b=',latex(f_b))
# Partial derivatives for Function b
partial_f_b_x1 = diff(f_b, x1)
partial_f_b_x2 = diff(f_b, x2)
partial_f_b_x3 = diff(f_b, x3)
print('dF_b1=',latex(partial_f_b_x1))
print('dF_b2=',latex(partial_f_b_x2))
print('dF_b3=',latex(partial_f_b_x3))
# Stationary points for Function b
stationary_points_b = solve((partial_f_b_x1, partial_f_b_x2, partial_f_b_x3), (x1, x2, x3))
print('sta_b=',latex(stationary_points_b))
# Hessian matrix for Function b
second_partial_f_b_x1x1 = diff(partial_f_b_x1, x1)
second_partial_f_b_x1x2 = diff(partial_f_b_x1, x2)
second_partial_f_b_x1x3 = diff(partial_f_b_x1, x3)
second_partial_f_b_x2x2 = diff(partial_f_b_x2, x2)
second_partial_f_b_x2x3 = diff(partial_f_b_x2, x3)
second_partial_f_b_x3x3 = diff(partial_f_b_x3, x3)
Hessian_b = Matrix([
    [second_partial_f_b_x1x1, second_partial_f_b_x1x2, second_partial_f_b_x1x3],
    [second_partial_f_b_x1x2, second_partial_f_b_x2x2, second_partial_f_b_x2x3],
    [second_partial_f_b_x1x3, second_partial_f_b_x2x3, second_partial_f_b_x3x3]
])
print('H_b=',latex(Hessian_b))
stationary_points_a, Hessian_a, stationary_points_b, Hessian_b

from sympy import N

# Function to check definiteness of Hessian matrix
def check_definiteness(matrix, variables, point):
    # Substitute the point into the matrix
    substituted_matrix = matrix.subs(zip(variables, point))

    # Evaluate eigenvalues
    eigenvalues = substituted_matrix.eigenvals()
    
    # All eigenvalues must be real for definiteness
    if not all([ev.is_real for ev in eigenvalues]):
        return "Indeterminate (complex eigenvalues)"

    # Check if eigenvalues are positive or negative
    positive, negative = False, False
    for ev in eigenvalues:
        if N(ev) > 0:
            positive = True
        elif N(ev) < 0:
            negative = True

    # Determine definiteness
    if positive and not negative:
        return "Positive definite (local minimum)"
    elif not positive and negative:
        return "Negative definite (local maximum)"
    elif positive and negative:
        return "Indefinite (saddle point)"
    else:
        return "Indeterminate"

# Check definiteness for each stationary point in Function a
definiteness_results_a = [(point, check_definiteness(Hessian_a, (x1, x2), point)) 
                          for point in stationary_points_a]

# Check definiteness for the stationary point in Function b
definiteness_result_b = check_definiteness(Hessian_b, (x1, x2, x3), (0, 0, 0))

#print(latex(definiteness_results_a))
print(latex(definiteness_result_b))


