import numpy as np
import sympy as sp
from sympy import latex

# Define symbolic variables
x1, x2 = sp.symbols('x1 x2')

# Define the objective function
F = x1**2 + x2**2 - 0.3*sp.cos(3*sp.pi*x1) - 0.4*sp.cos(4*sp.pi*x2) + 0.7

# Convert the symbolic expression to a Python function for numerical evaluation
F_num = sp.lambdify((x1, x2), F, 'numpy')

def simulated_annealing(F, bounds, initial_temp, cooling_rate, stopping_temp, max_iterations):
    """
    Simulated Annealing optimization for minimizing a given objective function.

    Parameters:
    - F: The objective function to minimize.
    - bounds: Tuple of (min, max) for x1 and x2.
    - initial_temp: Initial temperature for the annealing process.
    - cooling_rate: Rate at which the temperature will be decreased.
    - stopping_temp: Temperature at which the annealing process stops.
    - max_iterations: Maximum number of iterations per temperature level.

    Returns:
    - Best solution found and its objective function value.
    """
    # Random initial solution within the bounds
    current_solution = np.random.uniform(bounds[0], bounds[1], 2)
    current_value = F(*current_solution)
    best_solution = np.copy(current_solution)
    best_value = current_value
    
    temp = initial_temp

    # Main loop
    while temp > stopping_temp:
        for _ in range(max_iterations):
            # Generate a neighboring solution
            neighbor = current_solution + np.random.uniform(-0.1, 0.1, 2)
            # Ensure the neighbor is within bounds
            neighbor = np.clip(neighbor, bounds[0], bounds[1])
            neighbor_value = F(*neighbor)

            # Acceptance criterion (Metropolis criterion)
            if neighbor_value < current_value or np.random.rand() < np.exp((current_value - neighbor_value) / temp):
                current_solution = neighbor
                current_value = neighbor_value
                
                # Update best solution found
                if current_value < best_value:
                    best_solution = current_solution
                    best_value = current_value

        # Cooling down
        temp *= cooling_rate   #temp = temp*cooling_rate

    return best_solution, best_value

# Example of parameters setup for the simulated annealing
bounds = (-1, 1)
initial_temp = 10000
cooling_rate = 0.95
stopping_temp = 0.001
max_iterations = 1000

# Note: The function call is commented out to adhere to the instructions.
best_solution, best_value = simulated_annealing(F_num, bounds, initial_temp, cooling_rate, stopping_temp, max_iterations)
print("Best solution:", latex(best_solution))
print("Best objective function value:", latex(best_value))
