from sympy import symbols, diff, sqrt, latex, Matrix

# Define symbols
x1, x2 = symbols('x1 x2')

# Objective function
F = 4*(x1 - 5)**2 + (x2 - 6)**2

# Calculate gradient of F
gradient_F = [diff(F, x) for x in (x1, x2)]

print('F=',latex(F))
print('gF=',latex(gradient_F))
# Fletcher-Reeves method implementation
def fletcher_reeves(gradient, initial_guess, iterations=10):
    x = Matrix(initial_guess)  # Current guess
    grad_val = Matrix([g.subs(list(zip((x1, x2), x))) for g in gradient])  # Gradient at current guess
    p = -grad_val  # Initial search direction is the negative gradient
    for _ in range(iterations):
        # Perform line search (simplified here, normally requires more sophisticated approach)
        # For demonstration, let's assume alpha is a constant small step size for simplicity
        alpha = 0.01
        
        # Update guess
        x_new = x + alpha * p
        
        # Calculate new gradient
        grad_val_new = Matrix([g.subs(list(zip((x1, x2), x_new))) for g in gradient])
        
        # Calculate beta using Fletcher-Reeves formula
        beta = (grad_val_new.norm()**2) / (grad_val.norm()**2)
        
        # Update search direction
        p = -grad_val_new + beta * p
        
        # Update current guess and gradient
        x = x_new
        grad_val = grad_val_new
    
    return x

# Initial guess
initial_guess = (0, 0)

# Apply Fletcher-Reeves method
result = fletcher_reeves(gradient_F, initial_guess)

print('R',latex(result))
