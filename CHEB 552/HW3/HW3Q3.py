import numpy as np

def calculate_dependent_variables(x1, x7, x8):
    x4 = x1 * (1.12 + 0.13167 * x8 - 0.006667 * x8 ** 2)
    x5 = 1.22 * x4 - x1
    x2 = x1 * x8 - x5
    x6 = 89 + (x7 - (86.35 + 1.098 * x8 - 0.038 * x8 ** 2)) / 0.325
    x10 = -133 + 3 * x7
    x9 = 35.82 - 0.222 * x10
    x3 = 0.001 * x4 * x6 * x9 / (98 - x6)
    return x4, x5, x2, x6, x10, x3, x9

def satisfies_constraints(x4, x5, x2, x6, x10, x3, x9):
    if not (0 <= x3 <= 120): return False
    if not (0 <= x4 <= 5000): return False
    if not (0 <= x5 <= 2000): return False
    if not (0 <= x2 <= 16000): return False
    if not (85 <= x6 <= 93): return False
    if not (145 <= x10 <= 162): return False
    if not (1.2 <= x9 <= 4): return False
    return True

def objective_function(x1, x7, x8):
    x4, x5, x2, x6, x10, x3, x9 = calculate_dependent_variables(x1, x7, x8)
    if not satisfies_constraints(x4, x5, x2, x6, x10, x3, x9):
        return np.inf
    return -(0.063 * x4 * x7 - 5.04 * x1 - 0.035 * x2 - 10 * x3 - 3.36 * x5)

def luus_jaakola(x_bounds, max_iterations=1000, V_initial=0.1, V_reduction_factor=0.95):
    n = len(x_bounds)
    x = np.array([np.random.uniform(low, high) for low, high in x_bounds])
    best_value = objective_function(*x)
    best_solution = x
    
    for iteration in range(max_iterations):
        V = V_initial * (V_reduction_factor ** iteration)
        trial_x = x + np.random.uniform(-V, V, n)
        trial_x = np.clip(trial_x, [b[0] for b in x_bounds], [b[1] for b in x_bounds])
        
        trial_value = objective_function(*trial_x)
        if trial_value < best_value:
            best_value = trial_value
            best_solution = trial_x
            x = trial_x
    
    x4, x5, x2, x6, x10, x3, x9 = calculate_dependent_variables(*best_solution)
    dependent_vars = np.array([x2, x3, x4, x5, x6, x9, x10])
    return best_solution, best_value, dependent_vars

# Adjusted bounds for the independent variables (x1, x7, x8) based on the problem
x_bounds = [(1727, 1728), (94, 95), (10, 11)]

best_value = np.inf
for j in range(1000):
    solution, value, other_vars = luus_jaakola(x_bounds)
    if value != np.inf:
        best_value = value
        best_solution = solution
        best_other_vars = other_vars
        break

if best_value != np.inf:
    print("Best Solution: x1 = {:.5f}, x7 = {:.5f}, x8 = {:.5f}".format(*best_solution))
    print("Dependent Variables:", best_other_vars)
    print("Best Objective Function Value (Negative Profit): {:.5f}".format(best_value))
else:
    print("No feasible solution found within the given iterations.")
