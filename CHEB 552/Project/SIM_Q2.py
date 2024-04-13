from scipy.optimize import brentq
import numpy as np
from scipy.optimize import least_squares
HOAC = 3.2
O2 = 16
W_F_A0_data = np.array([6, 8.6, 3.3, 4, 4.9, 6.4, 9.6])
X_data = np.array([0.96, 0.96, 0.85, 0.91, 0.95, 0.97, 0.95]) 

def inverse_model_function_bounded(W_F_A0, k1, k2, k3):
    """
    Numerically solve for X given W_F_A0 using a bounded solver.
    """
    # Define the equation to solve
    def equation(X):
        return ((1 + k3 * O2)**2) * (k2 * HOAC * X - np.log(1 - X)) / (k1 * HOAC * O2) - W_F_A0
    
    # Use brentq for solving within bounds, providing more stability near edges
    try:
        X_solution = brentq(equation, 0.0001, 0.9999)  # Avoiding exact 0 or 1 to prevent log issues
        return X_solution
    except ValueError as e:
        print(f"No solution found for W_F_A0={W_F_A0} within the bounds: {e}")
        return np.nan  # Return NaN if no solution is found

def compute_X_for_all_data(W_F_A0_data, k1, k2, k3):
    """
    Compute X for all data points using the bounded model function.
    """
    computed_X = [inverse_model_function_bounded(w, k1, k2, k3) for w in W_F_A0_data]
    return computed_X

# Example usage with previously found best parameters
#params_example = [2.42877233, -0.24898018, 1.11437486]  # Use previously found best parameters


#computed_X_bounded = compute_X_for_all_data(W_F_A0_data, *params_example)
#print("Computed X values:", computed_X_bounded)
def objective_function_for_inverse_model(params, W_F_A0_data, X_data):
    """
    Objective function to minimize: the sum of squared differences between
    actual X data and computed X values from the inverse model.
    """
    k1, k2, k3 = params
    computed_X = np.array([inverse_model_function_bounded(w, k1, k2, k3) for w in W_F_A0_data])
    return np.sum((X_data - computed_X) ** 2)

def nelder_mead_simplex(obj_func, initial_simplex, tol=1e-6, max_iter=1000):
    """
    Nelder-Mead Simplex Algorithm for minimizing a scalar function.

    Parameters:
    - obj_func: The objective function to be minimized.
    - initial_simplex: Initial simplex given as a list of vectors (ndarray). There should be n + 1 vectors for an n-dimensional space.
    - tol: Tolerance for convergence.
    - max_iter: Maximum number of iterations.

    Returns:
    - Best parameter set found.
    """
    import numpy as np
    
    # Number of dimensions
    n = len(initial_simplex) - 1
    
    # Reflection, expansion, and contraction coefficients
    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5
    
    for iteration in range(max_iter):
        # Sort simplex by the objective function values
        initial_simplex.sort(key=obj_func)
        best = initial_simplex[0]
        worst = initial_simplex[-1]
        second_worst = initial_simplex[-2]

        # Compute the centroid of the simplex excluding the worst point
        centroid = np.sum(initial_simplex[:-1], axis=0) / n

        # Reflection
        reflected = centroid + alpha * (centroid - worst)
        if obj_func(best) <= obj_func(reflected) < obj_func(second_worst):
            initial_simplex[-1] = reflected
        elif obj_func(reflected) < obj_func(best):
            # Expansion
            expanded = centroid + gamma * (reflected - centroid)
            if obj_func(expanded) < obj_func(reflected):
                initial_simplex[-1] = expanded
            else:
                initial_simplex[-1] = reflected
        else:
            # Contraction
            contracted = centroid + rho * (worst - centroid)
            if obj_func(contracted) < obj_func(worst):
                initial_simplex[-1] = contracted
            else:
                # Shrink
                initial_simplex = [best + sigma * (x - best) for x in initial_simplex]

        # Check for convergence
        if np.std([obj_func(x) for x in initial_simplex]) < tol:
            break

    return best
def nelder_mead_optimize_for_params(W_F_A0_data, X_data, initial_guesses):
    """
    Use the Nelder-Mead algorithm to estimate the parameters k1, k2, k3.
    """
    # Creating an initial simplex
    initial_simplex = [np.array(initial_guesses)]
    for i in range(len(initial_guesses)):
        new_guess = np.array(initial_guesses)
        new_guess[i] += 0.1  # Arbitrary small variation to create initial simplex
        initial_simplex.append(new_guess)
    
    # Optimize using the Nelder-Mead Simplex algorithm
    result_params = nelder_mead_simplex(
        lambda params: objective_function_for_inverse_model(params, W_F_A0_data, X_data), 
        initial_simplex
    )
    return result_params

# Initial guesses for k1, k2, k3
initial_guesses = [1.21, 0.88, 0.36]  # Starting values for the parameters

# Running the optimization
#best_parameters = nelder_mead_optimize_for_params(W_F_A0_data, X_data, initial_guesses)
#best_parameters

def safe_log(x):
    """ Returns a safe logarithmic value to avoid math domain error. """
    if x <= 0:
        return -np.inf
    return np.log(x)

def refined_model_function(W_F_A0, k1, k2, k3):
    """ Adjusted model function for numerical stability. """
    def equation(X):
        if X <= 0 or X >= 1:
            return np.inf
        log_term = safe_log(1 - X)
        return ((1 + k3 * O2) ** 2) * (k2 * HOAC * X - log_term) / (k1 * HOAC * O2) - W_F_A0
    
    # Solve the equation using Brent's method within a safe range
    try:
        X_solution = brentq(equation, 0.01, 0.99)
        return X_solution
    except ValueError:
        return np.nan
def residuals(params, W_F_A0, X_data):
    k1, k2, k3 = params
    X_predicted = [refined_model_function(w, k1, k2, k3) for w in W_F_A0]
    return X_data - X_predicted

# Example use
initial_guesses = [1.0, 1.0, 1.0]
result = least_squares(residuals, initial_guesses, args=(W_F_A0_data, X_data))
print(result.x)  # Optimized parameters
import matplotlib.pyplot as plt

# Assume we have obtained best parameters, let's use the ones we obtained earlier
best_k1, best_k2, best_k3 = result.x[0],result.x[1],result.x[2]

# Calculate the estimated X values using the refined model function
estimated_X = [refined_model_function(w, best_k1, best_k2, best_k3) for w in W_F_A0_data]

# Plotting the experimental data and the estimated output
sorted_indices = np.argsort(W_F_A0_data)
sorted_W_F_A0 = W_F_A0_data[sorted_indices]
sorted_estimated_X = np.array(estimated_X)[sorted_indices]

# Re-plotting
plt.figure(figsize=(8, 5))
plt.scatter(W_F_A0_data, X_data, color='blue', label='Experimental X Data')
plt.plot(sorted_W_F_A0, sorted_estimated_X, 'r-', label='Estimated X')
plt.xlabel('W_F_A0')
plt.ylabel('X')
plt.title('Comparison of Experimental and Estimated X Values (Sorted)')
plt.legend()
plt.grid(True)
plt.show()
