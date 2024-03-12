import numpy as np

# Adjusted c_j values
c = 3.932 + np.array([-10.021, -21.096, -37.986, -9.846, -28.653, -18.918, -28.032, -14.640, -30.594, -26.111])

def objective_function(x_vars, epsilon=1e-6):
    x_2, x_5, x_6, x_7, x_8, x_9, x_10 = x_vars
    
    # Constraint
    x_4 = 1 - 2*x_5 - x_6 - x_7
    x_3 = 1 - x_7 - x_8 - 2*x_9 - x_10
    x_1 = 2 - 2*x_2 - 2*x_3 - x_6 - x_10
    x = np.array([x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9, x_10])
    XT = np.sum(x)
    # Ensure positive arguments for the logarithm
    log_args = np.where(x / XT + epsilon > 0, x / XT + epsilon, np.finfo(float).tiny)
    return np.sum(x * (c + np.log(log_args)))

def luus_jaakola(x_bounds, max_iterations=10000, V_initial=0.1, V_reduction_factor=0.999, epsilon=1e-6):
    n = len(x_bounds)
    x_vars = np.array([np.random.uniform(low, high) for low, high in x_bounds])  # Independent variables
    V = V_initial  # Initial search area

    for iteration in range(max_iterations):
        x_trial = x_vars + np.random.uniform(-V, V, n)
        x_trial = np.clip(x_trial, a_min=[low for low, _ in x_bounds], a_max=[high for _, high in x_bounds])
        
        if objective_function(x_trial, epsilon) < objective_function(x_vars, epsilon):
            x_vars = x_trial
        
        V *= V_reduction_factor  # Reduce search area
    x_2, x_5, x_6, x_7, x_8, x_9, x_10 = x_vars
    x_1 = 2 - 2*x_2 - 2*(1 - x_7 - x_8 - 2*x_9 - x_10) - x_6 - x_10
    x_3 = 1 - x_7 - x_8 - 2*x_9 - x_10
    x_4 = 1 - 2*x_5 - x_6 - x_7
    return np.array([x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9, x_10]), objective_function(x_vars, epsilon)

# Bounds for the independent variables (x2, x5, x6, x7, x8, x9, x10)
x_bounds = [(0, 1)] * 7  # Adjusted bounds to ensure positive values

# Execute the Luus-Jaakola method with initial xi = 0.2 and ri = 0.2 for i = 2, 5, 6, ..., 10
solution, value = luus_jaakola(x_bounds)
print("Solution:",  "{:.5f}".format(solution))
# print("Solution:", x_1,x_3,x_4)
print("Objective Function Value:", "{:.5f}".format(value))
