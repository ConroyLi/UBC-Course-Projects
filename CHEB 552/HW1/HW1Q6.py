from sympy import symbols, diff, Matrix, latex

# Define symbols
x1, x2 = symbols('x1 x2')

# Define the function f(x)
f = 100*(x2 - x1**2)**2 + (1 - x1)**2
print('f',latex(f))


# Compute the gradient vector
gradient = Matrix([diff(f, x1), diff(f, x2)])
print('gradient',latex(gradient))

# Compute the Hessian matrix
hessian = Matrix([[diff(gradient[0], x1), diff(gradient[0], x2)],
                  [diff(gradient[1], x1), diff(gradient[1], x2)]])
print('Hessian',latex(hessian))


# Evaluate the gradient at x* = [1, 1]
grad_at_x_star = gradient.subs({x1: 1, x2: 1})
print('grad_at_x_star',latex(grad_at_x_star))

# Check if x* = [1, 1] is a strong local minimum
# For this, the Hessian matrix at x* should be positive definite
hessian_at_x_star = hessian.subs({x1: 1, x2: 1})
print('hessian_at_x_star',latex(hessian_at_x_star))

# Check eigenvalues of the Hessian matrix at x* to determine if it's positive definite
eigenvalues_hessian_x_star = hessian_at_x_star.eigenvals()
print('eigenvalues_hessian_x_star',latex(eigenvalues_hessian_x_star))

eigenvalues = hessian.eigenvals()
print('eigenvalues',latex(eigenvalues))
gradient, hessian, grad_at_x_star, hessian_at_x_star, eigenvalues_hessian_x_star

