import numpy as np

# Example data (p, q) - user will replace this with actual data
#p_data = np.array([0, 0.3427, 0.5134, 0.6837, 0.8535, 1.0229, 1.1918, 1.5284, 1.7821, 2.1239, 2.5752, 3.2134, 3.5572, 3.679, 4.39, 4.8735, 5.1615])  # Input data
#q_data = np.array([0, 1.9267, 2.1711, 2.3449, 2.453, 2.5799, 2.6612, 2.7937, 2.851, 2.9159, 2.9991, 3.1062, 3.1438, 3.1551, 3.273, 3.3012, 3.3228])  # Actual output data
# p_data =np.array([0, 0.3559, 0.5332, 0.71, 0.8863, 1.0622, 1.2377, 1.4126, 1.5872, 1.8313, 2.1936, 2.7723, 3.2543, 3.6185, 3.759, 4.5065, 5.004, 5.2772])
# q_data =np.array([0, 1.2594, 1.4662, 1.6353, 1.7556, 1.8581, 1.9577, 2.0324, 2.1109, 2.2044, 2.3167, 2.4243, 2.4803, 2.5766, 2.6024, 2.7204, 2.7979, 2.813])
# p_data =np.array([0, 0.3651, 0.5470, 0.7284, 0.9093, 1.0898, 1.2697, 1.4493, 1.6283, 1.8663, 2.2458, 2.8411, 3.2951, 3.6798, 3.8230, 4.6036, 5.0910, 5.3930])
# q_data =np.array([0, 0.5874, 0.7237, 0.8365, 0.9305, 0.9803, 1.0818, 1.1316, 1.2030, 1.2768, 1.3440, 1.4685, 1.5268, 1.6034, 1.6198, 1.7199, 1.7857, 1.8120])
# p_data =np.array([0, 0.3849, 0.5766, 0.7678, 0.9586, 1.1488, 1.3385, 1.5277, 1.7165, 1.9576, 2.3329, 2.9556, 3.4313, 3.8025, 3.9350, 4.7008, 5.2216, 5.5087])
# q_data =np.array([0, 0.3102, 0.3853, 0.4525, 0.4944, 0.5503, 0.5954, 0.6297, 0.6767, 0.7265, 0.7810, 0.8623, 0.9253, 0.9662, 0.9826, 1.0526, 1.0982, 1.1194])
p_data =np.array([0, 0.3981, 0.5964, 0.7941, 0.9914, 1.1881, 1.3843, 1.5801, 1.7753, 2.0347, 2.4199, 3.0243, 3.5130, 3.8945, 4.0469, 4.7785, 5.2651, 5.5550])
q_data =np.array([0, 0.1175, 0.1551, 0.1833, 0.2054, 0.2162, 0.2547, 0.2693, 0.2923, 0.3163, 0.3407, 0.3839, 0.4145, 0.4370, 0.4474, 0.4737, 0.5066, 0.5202])

def algebraic_model(p, q_sat, k, t):
    return (q_sat * k * p) / ((1 + (k * p) ** t) ** (1/t))

def objective_function(params):
    q_sat, k, t = params
    predicted_q = algebraic_model(p_data, q_sat, k, t)
    return np.sum((q_data - predicted_q) ** 2)

# Initial guesses for q_sat, k, t
#initial_guesses = [4, 17, 0.6]
#initial_guesses = [6, 9, 0.5]
#initial_guesses = [5, 2, 0.4]
#initial_guesses = [5.5, 0.7, 0.3]
initial_guesses = [3, 0.3, 0.3] # Initial guesses for q_sat, k, t
parameter_variations = [0.5, 1, 0.04]

# Creating an initial simplex
# For n parameters, we need n + 1 initial guesses
adjusted_initial_simplex = [np.array(initial_guesses)]
for i in range(len(initial_guesses)):
    new_guess = np.array(initial_guesses)
    new_guess[i] += parameter_variations[i]  # Apply different variations based on parameter magnitude
    adjusted_initial_simplex.append(new_guess)

def nelder_mead_simplex(obj_func, initial_simplex, tol=1e-6, max_iter=1000):
   
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

best_parameters = nelder_mead_simplex(objective_function, adjusted_initial_simplex, tol=1e-6, max_iter=1000)
print("Estimated Parameters:", np.round(best_parameters,2))

from sklearn.utils import resample
import numpy as np

def bootstrap_CI(data, model_func, optimization_func, initial_guesses, n_iterations=1000, ci=95):

    np.random.seed(42)  # For reproducibility
    p_data, q_data = data
    estimates = []

    for _ in range(n_iterations):
        # Bootstrap sample
        bs_indices = np.random.choice(range(len(p_data)), size=len(p_data), replace=True)
        bs_p_data = p_data[bs_indices]
        bs_q_data = q_data[bs_indices]
        
        # Update objective function for bootstrap sample
        def bs_objective_function(params):
            predicted_q = model_func(bs_p_data, *params)
            return np.sum((bs_q_data - predicted_q) ** 2)
        
        # Estimate parameters for bootstrap sample
        bs_estimates = optimization_func(bs_objective_function, initial_guesses)
        estimates.append(bs_estimates)
    
    # Calculate confidence intervals
    estimates = np.array(estimates)
    lower_bounds = np.percentile(estimates, (100 - ci) / 2, axis=0)
    upper_bounds = np.percentile(estimates, 100 - (100 - ci) / 2, axis=0)
    
    return list(zip(lower_bounds, upper_bounds))

# Note: The actual execution of bootstrap_CI is commented out to prevent execution.
# The optimization_func should be adapted to return the best parameters from the Nelder-Mead implementation.
ci_results = bootstrap_CI((p_data, q_data), algebraic_model, nelder_mead_simplex, adjusted_initial_simplex, n_iterations=1000, ci=95)
print("95% Confidence Intervals:", np.round(ci_results,2))
