from sympy import symbols, diff, hessian, Matrix, latex, solve
from sympy import init_printing

# Part a

# Initialize printing for nice symbolic representation
init_printing(use_unicode=True)

# Define symbols
x1, x2, x3, x4 = symbols('x1 x2 x3 x4')

# Function 1
F1 = 1 + x1 + x2 + x3 + x4 + x1*x2 + x1*x3 + x1*x4 + x2*x3 + x2*x4 + x3*x4 + x1**2 + x2**2 + x3**2 + x4**2

# Calculate gradient of F1
gradient_F1 = [diff(F1, x) for x in (x1, x2, x3, x4)]

# Calculate Hessian matrix of F1
Hessian_F1 = hessian(F1, (x1, x2, x3, x4))

# Function 2
x1, x2 = symbols('x1 x2')  # Redefine symbols for function 2
F2 = 8*x1**2 + 4*x1*x2 + 5*x2**2

# Calculate gradient of F2
gradient_F2 = [diff(F2, x) for x in (x1, x2)]

# Calculate Hessian matrix of F2
Hessian_F2 = hessian(F2, (x1, x2))

print('f1=',latex(F1))
print('H1=',latex(Hessian_F1))
print('f2=',latex(F2))
print('H2=',latex(Hessian_F2))
# Newton's method implementation
def newtons_method(gradient, hessian, initial_guess, iterations=10):
    x = Matrix(initial_guess)
    for _ in range(iterations):
        grad_val = Matrix([g.subs(list(zip((x1, x2, x3, x4), x))) for g in gradient])
        hess_val = hessian.subs(list(zip((x1, x2, x3, x4), x)))
        x -= hess_val.inv() * grad_val
    return x

# Initial guesses for Function 1
initial_guesses_F1 = [(-3, -30, -4, -0.1), (0.5, 1.0, 8.0, -0.7)]

# Apply Newton's method for Function 1 with both initial guesses
results_F1 = [newtons_method(gradient_F1, Hessian_F1, guess) for guess in initial_guesses_F1]

# Initial guess for Function 2
initial_guess_F2 = (10, 10)

# Apply Newton's method for Function 2
result_F2 = newtons_method(gradient_F2, Hessian_F2, initial_guess_F2, iterations=10)

print('R1=', latex(results_F1))
print('R2=', latex(result_F2))