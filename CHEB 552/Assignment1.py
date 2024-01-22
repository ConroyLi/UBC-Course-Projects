from sympy import symbols, diff, solve, hessian, Matrix

# Define symbols
x1, x2, x = symbols('x1 x2 x')

# Function 1
F_a = (x1 - 2)**4 + (x1 - 2*x2)**2

# Partial derivatives for Function a
partial_F_a_x1 = diff(F_a, x1)
partial_F_a_x2 = diff(F_a, x2)

# Stationary points for Function a
stationary_points_a = solve((partial_F_a_x1, partial_F_a_x2), (x1, x2))

# Function b
F_b = 2*x1**3 + x2**2 + x1**2*x2**2 + 4*x1*x2 + 3

# Partial derivatives for Function b
partial_F_b_x1 = diff(F_b, x1)
partial_F_b_x2 = diff(F_b, x2)

# Stationary points for Function b
stationary_points_b = solve((partial_F_b_x1, partial_F_b_x2), (x1, x2))

# Function c
F_c = 12*x**5 - 45*x**4 + 40*x**3 + 5

# Derivative for Function c
F_c_prime = diff(F_c, x)

# Stationary points for Function c
stationary_points_c = solve(F_c_prime, x)

print(stationary_points_a, stationary_points_b, stationary_points_c)


# Second derivatives for Function a
second_partial_F_a_x1x1 = diff(partial_F_a_x1, x1)
second_partial_F_a_x1x2 = diff(partial_F_a_x1, x2)
second_partial_F_a_x2x2 = diff(partial_F_a_x2, x2)

# Hessian matrix for Function a
Hessian_a = Matrix([[second_partial_F_a_x1x1, second_partial_F_a_x1x2], 
                    [second_partial_F_a_x1x2, second_partial_F_a_x2x2]])

# Evaluate the Hessian matrix at the stationary point for Function a
Hessian_a_at_stationary_a = Hessian_a.subs({x1: 2, x2: 1})

# Second derivatives for Function b
second_partial_F_b_x1x1 = diff(partial_F_b_x1, x1)
second_partial_F_b_x1x2 = diff(partial_F_b_x1, x2)
second_partial_F_b_x2x2 = diff(partial_F_b_x2, x2)

# Hessian matrix for Function b
Hessian_b = Matrix([[second_partial_F_b_x1x1, second_partial_F_b_x1x2], 
                    [second_partial_F_b_x1x2, second_partial_F_b_x2x2]])

# Evaluate the Hessian matrix at the stationary point for Function b
Hessian_b_at_stationary_b = Hessian_b.subs({x1: 0, x2: 0})

# Second derivative for Function c
F_c_double_prime = diff(F_c_prime, x)

# Evaluate the second derivative at the stationary points for Function c
second_derivative_at_stationary_c = [F_c_double_prime.subs(x, point) for point in stationary_points_c]

print(Hessian_a_at_stationary_a, Hessian_b_at_stationary_b, second_derivative_at_stationary_c)

# Checking the definiteness of the Hessian matrices

def check_definiteness(matrix):
    """ Check the definiteness of a matrix """
    eigenvalues = matrix.eigenvals()
    if all(value > 0 for value in eigenvalues):
        return "Positive definite (local minimum)"
    elif all(value < 0 for value in eigenvalues):
        return "Negative definite (local maximum)"
    else:
        return "Indefinite (saddle point)"

# Definiteness of the Hessian matrix for Function a and b
definiteness_a = check_definiteness(Hessian_a_at_stationary_a)
definiteness_b = check_definiteness(Hessian_b_at_stationary_b)

definiteness_a, definiteness_b

