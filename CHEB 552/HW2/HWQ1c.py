from sympy import symbols, diff, sqrt, latex, Matrix
# Define symbols
x1, x2 = symbols('x1 x2')

# Objective function
F = x1 - x2 + 2*x1**2 + 2*x1*x2 + x2**2

# Calculate gradient of F
gradient_F = [diff(F, x) for x in (x1, x2)]
print('F=',latex(F))
print('gF=',latex(gradient_F))
# DFP method implementation
def dfp(gradient, initial_guess, iterations=10):
    x = Matrix(initial_guess)  # Current guess
    H = Matrix.eye(2)  # Initial Hessian approximation (identity matrix)
    for _ in range(iterations):
        grad_val = Matrix([g.subs(list(zip((x1, x2), x))) for g in gradient])  # Gradient at current guess
        p = -H * grad_val  # Search direction
        
        # For simplicity, assume a small constant step size
        alpha = 0.1
        
        # Update guess
        x_new = x + alpha * p
        
        # Calculate new gradient
        grad_val_new = Matrix([g.subs(list(zip((x1, x2), x_new))) for g in gradient])
        
        # Calculate differences
        delta_x = x_new - x
        delta_grad = grad_val_new - grad_val
        
        # Update Hessian approximation using DFP formula
        dx_dxT = delta_x * delta_x.transpose()
        Hdg_dgT = H * delta_grad * delta_grad.transpose() * H
        dxT_dg = delta_x.transpose() * delta_grad
        HdgT_dg = delta_grad.transpose() * H * delta_grad
        
        H += dx_dxT / dxT_dg[0] - Hdg_dgT / HdgT_dg[0]
        
        # Update current guess and gradient
        x = x_new
        grad_val = grad_val_new
    
    return x

# Initial guess
initial_guess = (0, 0)

# Apply DFP method
result = dfp(gradient_F, initial_guess)

print('R',latex(result))