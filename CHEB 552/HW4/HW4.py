import scipy.stats
import numpy as np
import matplotlib.pyplot as plt

#GNM Method
def objective_function(params, X, Y, model):
   
    predictions = model(X, params)
    residuals = Y - predictions
    return np.sum(residuals**2)

def jacobian_matrix(params, X, Y, model, epsilon=1e-6):
   
    jacobian = np.zeros((len(Y), len(params)))
    for i in range(len(params)):
        params_eps = params.copy()
        params_eps[i] += epsilon
        grad_approx = (model(X, params_eps) - model(X, params)) / epsilon
        jacobian[:, i] = grad_approx
    return jacobian

def gauss_newton_marquardt(X, Y, model, initial_params, alpha=0.05, max_iterations=1000, tau=9e-4, epsilon=1e-6, convergence_threshold=1e-4):
    
    params = initial_params.copy()
    iteration_details = []

    n = len(Y)
    p = len(initial_params)
    for iteration in range(max_iterations):
        J = jacobian_matrix(params, X, Y, model, epsilon)
        residuals = Y - model(X, params)
        obj_function_value = np.sum(residuals**2)
        H = J.T @ J + tau * np.eye(len(params))  # Hessian approximation with Marquardt's modification
        g = J.T @ residuals
        param_update = np.linalg.inv(H) @ g
        params += param_update

        iteration_details.append({
            'iteration': iteration + 1,
            'objective_function_value': obj_function_value,
            'parameters': params.copy()
        })

        if np.linalg.norm(param_update) < convergence_threshold:
            print(f"Convergence reached after {iteration+1} iterations.")
            break
    
    # Calculate standard errors and confidence intervals
    sigma_squared_hat = np.sum(residuals**2) / (n - p)
    parameter_variances = np.diag(sigma_squared_hat * np.linalg.inv(H))
    standard_errors = np.sqrt(parameter_variances)
    t_value = scipy.stats.t.ppf(1 - alpha / 2, n - p)
    confidence_intervals = np.array([params - t_value * standard_errors, params + t_value * standard_errors]).T
    
    return params, confidence_intervals, iteration_details

# Question 1

# Model A function definition
def model_A(X, params):
    
    kH, kR, KA = params
    pA = X  # Assuming X is just pA for simplification
    
    R = kH + ((kH**2) / (2 * kR)) * ((1 + KA * pA)**2) / (KA * pA)
    rAi = R - np.sqrt(R**2 - kH**2)
    return rAi

# Data preparation (simplified example)
pA_600 = np.array([1.0, 7.0, 4.0, 10.0, 14.6, 5.5, 8.5, 3.0, 0.22, 1.0])
pA_575 = np.array([1.0, 3.0, 5.0, 7.0, 9.6])
rAi_600 = np.array([0.0392, 0.0416, 0.0416, 0.0326, 0.0247, 0.0415, 0.0376, 0.0420, 0.0295, 0.0410])
rAi_575 = np.array([0.0227, 0.0277, 0.0255, 0.0217, 0.0183]) 

# Initial parameters from Table 1 for Model A at 600 oF (for example)
initial_params_A_600 = np.array([7e-2, 70e-2, 50e-2])  # Assuming kH, kR, KA are given in correct units
initial_params_A_575 = np.array([7e-2, 20e-2, 40e-2])
# Using the gauss_newton_marquardt function to estimate parameters for Model A
estimated_params_A_600, CI_A_600, iteration_details_A_600 = gauss_newton_marquardt(pA_600, rAi_600, model_A, initial_params_A_600)
estimated_params_A_575, CI_A_575, iteration_details_A_575 = gauss_newton_marquardt(pA_575, rAi_575, model_A, initial_params_A_575)

estimated_rAi_600 = model_A(pA_600, estimated_params_A_600)
estimated_rAi_575 = model_A(pA_575, estimated_params_A_575)

# Model B function definition
def model_B(X, params):
    
    kH, kR, KA = params
    lambda_ = -0.7
    pA = X  # Assuming X is just pA for simplification
    
    rAi = ((kH**lambda_) + (kR * KA * pA / ((1 + KA * pA)**2))**lambda_)**(1/ lambda_)
    return rAi

# Initial parameters from Table 2 for Model B at 600 oF (as an example)
initial_params_B = np.array([9e-2, 60e-2, 50e-2]) 

# Using the gauss_newton_marquardt function to estimate parameters for Model B
estimated_params_B_600, CI_B_600, iteration_details_B_600 = gauss_newton_marquardt(pA_600, rAi_600, model_B, initial_params_B)
estimated_params_B_575, CI_B_575, iteration_details_B_575 = gauss_newton_marquardt(pA_575, rAi_575, model_B, initial_params_B)

estimated_rBi_600 = model_A(pA_600, estimated_params_B_600)
estimated_rBi_575 = model_A(pA_575, estimated_params_B_575)

# Postprocessing
'''
print('A600=', estimated_params_A_600, CI_A_600)
print('A575=', estimated_params_A_575, CI_A_575)
print('B600=',estimated_params_B_600, CI_B_600)
print('B575=',estimated_params_B_575, CI_B_575)
'''
'''
for detail in iteration_details_B_575:
    iteration_row = f"{detail['iteration']} & {detail['objective_function_value']:.6f} & "
    iteration_row += " & ".join([f"{param:.3f}" for param in detail['parameters']])
    iteration_row += " \\\\"
    print(iteration_row)
'''

# Plotting the original vs. predicted rAi values
plt.figure(figsize=(20, 10))
plt.subplot(2, 2, 1)
plt.plot(pA_600, rAi_600, 'o', label='Original rAi values for T = 600', markersize=8)
plt.plot(pA_600, estimated_rAi_600, 'x', label='Predicted rAi values with Estimated Parametersfor T = 600', markersize=8)
plt.xlabel('Partial Pressure of sec-butyl alcohol (pA) for T = 600', fontsize=14)
plt.ylabel('Initial Rate (rAi) for T = 600', fontsize=14)
plt.title('Comparison of Original and Predicted Initial Rates (rAi)  for T = 600', fontsize=16)
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 2)
plt.plot(pA_575, rAi_575, 'o', label='Original rAi values for T = 575', markersize=8)
plt.plot(pA_575, estimated_rAi_575, 'x', label='Predicted rAi values with Estimated Parametersfor T = 575', markersize=8)
plt.xlabel('Partial Pressure of sec-butyl alcohol (pA) for T = 575', fontsize=14)
plt.ylabel('Initial Rate (rAi) for T = 575', fontsize=14)
plt.title('And for T = 575 for Model A', fontsize=16)
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 3)
plt.plot(pA_600, rAi_600, 'o', label='Original rAi values for T = 600', markersize=8)
plt.plot(pA_600, estimated_rBi_600, 'x', label='Predicted rAi values with Estimated Parametersfor T = 600', markersize=8)
plt.xlabel('Partial Pressure of sec-butyl alcohol (pA) for T = 600', fontsize=14)
plt.ylabel('Initial Rate (rAi) for T = 600', fontsize=14)
plt.title('Comparison of Original and Predicted Initial Rates (rAi)  for T = 600', fontsize=16)
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 4)
plt.plot(pA_575, rAi_575, 'o', label='Original rAi values for T = 575', markersize=8)
plt.plot(pA_575, estimated_rBi_575, 'x', label='Predicted rAi values with Estimated Parametersfor T = 575', markersize=8)
plt.xlabel('Partial Pressure of sec-butyl alcohol (pA) for T = 575', fontsize=14)
plt.ylabel('Initial Rate (rAi) for T = 575', fontsize=14)
plt.title('And for T = 575 for Model B', fontsize=16)
plt.legend()
plt.grid(True)

plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.4, hspace=0.2)
plt.tight_layout()
# plt.show()

# Question 2
# Define the model function for the Oxidation of Propylene problem
def oxidation_propylene_model(X, params):
    
    ko, kp = params
    co, cp, n = X.T  # Transpose to unpack the column-wise data into variables
    
    rp = (ko * kp * co**0.5 * cp) / (ko * co**0.5 + n * kp * cp)
    return rp

# Prepare the data based on the table provided for Problem 2
# Simplified example data - in practice, you should input the exact data from the PDF
co_values = np.array([3.07, 3.18, 1.24, 3.85, 3.15, 3.89, 6.48, 3.13, 3.14, 7.93, 7.79, 8.03, 7.78, 3.03, 8.00, 8.22, 6.13, 8.41, 7.75, 3.10, 1.25, 7.89, 3.06])  # Oxygen concentration
cp_values = np.array([3.05, 1.37, 3.17, 3.02, 4.31, 2.78, 3.11, 2.96, 2.84, 1.46, 1.38, 1.42, 1.49, 3.01, 1.35, 1.52, 5.95, 1.46, 5.68, 1.36, 1.42, 3.18, 2.87])  # Propylene concentration
n_values = np.array([0.658, 0.439, 0.452, 0.695, 0.635, 0.67, 0.76, 0.642, 0.665, 0.525, 0.483, 0.522, 0.53, 0.635, 0.48, 0.544, 0.893, 0.517, 0.996, 0.416, 0.367, 0.835, 0.609])  # Stoichiometric number
rp_values = np.array([2.73, 2.86, 3.00, 2.64, 2.60, 2.73, 2.56, 2.69, 2.77, 2.91, 2.87, 2.97, 2.93, 2.75, 2.9, 2.94, 2.38, 2.89, 2.41, 2.81, 2.86, 2.59, 2.76])  # Rate of propylene disappearance

# Combine co, cp, and n into a single array for the model input
X_data = np.column_stack((co_values, cp_values, n_values))

# Initial parameters based on the hints provided
initial_params_oxidation = np.array([1330, 0.5])  # Initial guesses for ko and kp

# Use the gauss_newton_marquardt function to estimate parameters for the Oxidation of Propylene
estimated_params_oxidation,CI_Q2,iteration_details_Q2 = gauss_newton_marquardt(X_data, rp_values, oxidation_propylene_model, initial_params_oxidation)
estimated_oxidation = oxidation_propylene_model(X_data, estimated_params_oxidation)
# print('Q2=',estimated_params_oxidation,CI_Q2)

plt.figure(figsize=(10, 15))
plt.subplot(3,1,1)
plt.plot(co_values, rp_values, 'o', label='Original oxidation values ', markersize=8)
plt.plot(co_values, estimated_oxidation, 'x', label='Predicted oxidation values with Estimated Parameters', markersize=8)
plt.xlabel('co_values', fontsize=14)
plt.ylabel('Oxidation values', fontsize=14)
plt.title('Comparison of Original and Predicted Oxidation values', fontsize=16)
plt.legend()
plt.grid(True)

plt.subplot(3,1,2)
plt.plot(cp_values, rp_values, 'o', label='Original oxidation values ', markersize=8)
plt.plot(cp_values, estimated_oxidation, 'x', label='Predicted oxidation values with Estimated Parameters', markersize=8)
plt.xlabel('cp_values', fontsize=14)
plt.ylabel('Oxidation values', fontsize=14)
plt.title('Comparison of Original and Predicted Oxidation values', fontsize=16)
plt.legend()
plt.grid(True)
plt.subplot(3,1,3)
plt.plot(n_values, rp_values, 'o', label='Original oxidation values ', markersize=8)
plt.plot(n_values, estimated_oxidation, 'x', label='Predicted oxidation values with Estimated Parameters', markersize=8)
plt.xlabel('n_values', fontsize=14)
plt.ylabel('Oxidation values', fontsize=14)
plt.title('Comparison of Original and Predicted Oxidation values', fontsize=16)
plt.legend()
plt.grid(True)


plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.4, hspace=0.2)
plt.tight_layout()
# plt.show()

plt.figure(figsize=(10, 8))
plt.plot(rp_values, estimated_oxidation, 'o', label='Original oxidation values ', markersize=8)
#plt.plot(co_values, estimated_oxidation, 'x', label='Predicted oxidation values with Estimated Parameters', markersize=8)
plt.xlabel('Original Data', fontsize=14)
plt.ylabel('Estimated Data', fontsize=14)
plt.title('Comparison of Original and Predicted Oxidation values', fontsize=16)
plt.legend()
plt.grid(True)
# plt.show()


# LJ 
def luus_jaakola_optimize(X, Y, model, param_bounds, initial_vicinity, reduction_factor, max_iterations):
    num_params = len(param_bounds)
    current_params = np.array([np.random.uniform(low, high) for low, high in param_bounds])
    best_obj = np.inf
    
    for iteration in range(max_iterations):
        vicinity = initial_vicinity * (reduction_factor ** iteration)
        sample_params = current_params + np.random.uniform(-vicinity, vicinity, num_params)
        
        # Ensure sample_params are within bounds
        sample_params = np.clip(sample_params, [low for low, _ in param_bounds], [high for _, high in param_bounds])
        
        obj = np.sum((Y - model(X, sample_params))**2)  # Objective function: sum of squared residuals
        
        if obj < best_obj:
            best_obj = obj
            current_params = sample_params
    
    return current_params

param_bounds_A = [(0.001, 0.2),  # Bounds for kH
                  (0.001, 1), # Bounds for kR
                  (0.001, 1)] # Bounds for KA
param_bounds_problem2 = [(1300, 1400),  # Bounds for ko
                          (0.5, 0.7)]    # Bounds for kp
# Initial vicinity - Start with a broad search range
initial_vicinity_A = 0.01

# Reduction factor - Reduce the vicinity size by this factor each iteration
reduction_factor_A = 0.95

# Maximum number of iterations
max_iterations_A = 10000

estimated_params_A_600_LJ = luus_jaakola_optimize(pA_600, rAi_600, model_A, param_bounds_A, initial_vicinity_A, reduction_factor_A, max_iterations_A)
estimated_params_A_575_LJ = luus_jaakola_optimize(pA_575, rAi_575, model_A, param_bounds_A, initial_vicinity_A, reduction_factor_A, max_iterations_A)
estimated_params_B_600_LJ = luus_jaakola_optimize(pA_600, rAi_600, model_B, param_bounds_A, initial_vicinity_A, reduction_factor_A, max_iterations_A)
estimated_params_B_575_LJ = luus_jaakola_optimize(pA_575, rAi_575, model_B, param_bounds_A, initial_vicinity_A, reduction_factor_A, max_iterations_A)
estimated_params_oxidation_LJ = luus_jaakola_optimize(X_data, rp_values, oxidation_propylene_model, param_bounds_problem2, initial_vicinity_A, reduction_factor_A, max_iterations_A)


print('A600=', estimated_params_A_600_LJ)
print('A575=', estimated_params_A_575_LJ)
print('B600=',estimated_params_B_600_LJ)
print('B575=',estimated_params_B_575_LJ)
print('Q2=',estimated_params_oxidation_LJ)