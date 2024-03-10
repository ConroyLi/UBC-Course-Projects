import numpy as np
import sympy as sp
from sympy import latex

# Define symbolic variables
x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')

def F1(x):
    """Objective function F(x) = 100(x2 - x1^2)^2 + (1 - x1)^2."""
    x1, x2 = x
    return 100 * (x2 - x1**2)**2 + (1 - x1)**2
def F2(x):
    """Objective function F(x) = 100(x2 - x1^2)^2 + (1 - x1)^2."""
    x1, x2, x3, x4 = x
    return (x1 + 10*x2)**2 + 5*(x3 - x4)**2 + (x2 - 2*x3)**4 + 10*(x1 - x4)**4

def nelder_mead(f, initial_point, alpha=1, gamma=2, rho=0.5, sigma=0.5, max_iterations=1000, tolerance=1e-10):
    """
    Nelder-Mead optimization algorithm.

    Parameters:
    - f: The objective function to minimize.
    - initial_point: Starting point for the optimization.
    - alpha, gamma, rho, sigma: Parameters for reflection, expansion, contraction, and shrink.
    - max_iterations: Maximum number of iterations before stopping.
    - tolerance: Convergence tolerance.

    Returns:
    - The minimum point found and the function value at that point.
    """
    # Step 1: Initialization
    n = len(initial_point)
    simplex = [initial_point]
    for i in range(n):
        x = initial_point.copy()
        x[i] += 1  # Small perturbation
        simplex.append(x)
    
    for iteration in range(max_iterations):
        simplex.sort(key=f)
        f_values = [f(x) for x in simplex]

        if np.max(f_values) - np.min(f_values) < tolerance:
            return simplex[0], f(simplex[0])
        
        centroid = np.sum(simplex[:-1], axis=0) / n
        reflected = centroid + alpha * (centroid - simplex[-1])
        f_reflected = f(reflected)

        if f_values[0] <= f_reflected < f_values[-2]:
            simplex[-1] = reflected
            continue

        if f_reflected < f_values[0]:
            expanded = centroid + gamma * (reflected - centroid)
            if f(expanded) < f_reflected:
                simplex[-1] = expanded
            else:
                simplex[-1] = reflected
            continue

        contracted = centroid + rho * (simplex[-1] - centroid)
        if f(contracted) < f_values[-1]:
            simplex[-1] = contracted
            continue

        # Contraction
        contracted = centroid + rho * (simplex[-1] - centroid)
        if f(contracted) < f_values[-1]:
            simplex[-1] = contracted
            continue

        # Shrink
        simplex = [simplex[0]] + [simplex[0] + sigma * (x - simplex[0]) for x in simplex[1:]]

    # If max_iterations were reached without convergence
    return simplex[0], f(simplex[0])

# Starting point
initial_point_1 = [-1.2, 1]
initial_point_2 = [-3, -1, 0, 1]

# Run the Nelder-Mead optimization
minimum_point_1, minimum_value_1 = nelder_mead(F1, initial_point_1)
minimum_point_2, minimum_value_2 = nelder_mead(F2, initial_point_2)
print("Minimum point:", minimum_point_1)
print("Minimum value:", minimum_value_1)
print("Minimum point:", minimum_point_2)
print("Minimum value:", minimum_value_2)
