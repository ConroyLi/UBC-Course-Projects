import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt
# Constants given in the problem
HOAC = 3.2
O2 = 16

# Initial parameter valu
initial_params = np.array([1, 1, 1], dtype=np.float64)
# np.array([10**4, 10**9, 10**6]
# Data provided by the user
# W_F_A0_data = np.array([6, 8.6, 3.3, 4, 4.9, 6.4, 9.6])
# X_data = np.array([0.96, 0.96, 0.85, 0.91, 0.95, 0.97, 0.95]) # Convert percentages to fractions for computation
W_F_A0_data = np.array([2.8, 3.3, 4, 5.3, 7.2])
X_data = np.array([0.37, 0.45, 0.51, 0.69, 0.66])
# Function to calculate the theoretical W/F_A0 based on the parameters and X
def model_function(X, k1, k2, k3):
    return ((1 + k3 * O2)**2 * (k2 * HOAC * X - np.log(1-X)) / (k1 * HOAC * O2))
#A = model_function(0.96, 1, 1, 1)
#print(A)

# Function to calculate the residuals (difference between actual data and model predictions)
def residuals(params, X_data, W_F_A0_data):
    k1, k2, k3 = params
    return W_F_A0_data - model_function(X_data, k1, k2, k3)

# Function to calculate the Jacobian matrix for the Gauss-Newton method
def jacobian(X, k1, k2, k3):
    # Partial derivatives of the model function with respect to each parameter
    J = np.empty((X.size, 3))
    J[:, 0] = -(O2*k3 + 1)**2*(X*HOAC*k2 - np.log(1 - X))/(HOAC*O2*k1**2)
    J[:, 1] = X*(O2*k3 + 1)**2/(O2*k1)
    J[:, 2] = 2*(O2*k3 + 1)*(X*HOAC*k2 - np.log(1 - X))/(HOAC*k1)
    return J

# Gauss-Newton method with Marquardt's modification
def gauss_newton_marquardt(X_data, W_F_A0_data, initial_params, max_iterations=100, lambda_factor=10, tol=1e-6):
    params = initial_params
    n = len(W_F_A0_data)
    lambda_param = np.identity(3)
    iteration_details = []
    for i in range(max_iterations):
        res = residuals(params, X_data, W_F_A0_data)
        J = jacobian(X_data, *params)
        
        # Calculate the parameter update with Marquardt's modification
        update = np.linalg.inv(J.T @ J + lambda_param * np.max(np.diagonal(J.T @ J)) * lambda_factor) @ J.T @ res
        
        # Check for convergence (if the update is very small, break)
        if np.linalg.norm(update) < tol:
            break
        
        # Update the parameters
        params += update
        
        # Output information required by the user
        objective_function_value = np.sum(res**2)
        iteration_details.append((i + 1, objective_function_value, params.copy()))
    # Calculating the covariance matrix and 95% confidence intervals for parameters
    J = jacobian(X_data, *params)  # Recompute the Jacobian with final parameters
    covariance_matrix = np.linalg.pinv(J.T @ J)
    alpha = 0.05  # 95% confidence level
    dof = n - len(params)  # degrees of freedom
    t_critical = scipy.stats.t.ppf(1 - alpha/2, dof)  # t-distribution critical value
    parameter_errors = t_critical * np.sqrt(np.diagonal(covariance_matrix) / n)
    confidence_intervals = [(params[i] - parameter_errors[i], params[i] + parameter_errors[i]) for i in range(len(params))]

    return params, confidence_intervals,iteration_details

# Uncomment the line below to run the Gauss-Newton method with the provided data
estimated_params, confidence_intervals,iteration_details = gauss_newton_marquardt(X_data, W_F_A0_data, initial_params)
print(np.round(estimated_params,2))
print(np.round(confidence_intervals,2))
# NOTE: The scipy.stats module is needed for the calculation of the t-distribution critical value,
# but since we don't have internet access, we cannot import it here. This part of the code will not execute.
# In an actual environment with internet access, we would uncomment the import statement at the beginning and
# include 'import scipy.stats' to make this part work.

# Assuming 'model_function' is defined and we have the estimated parameters 'estimated_params'
# Let's use the estimated parameters to calculate the estimated W/F_A0 values
estimated_W_F_A0 = model_function(X_data, *estimated_params)

# Plotting the estimated vs. original values of W/F_A0
plt.figure(figsize=(8, 6))
plt.scatter(W_F_A0_data, estimated_W_F_A0, color='blue', label='Estimated vs. Original')
plt.plot(W_F_A0_data, W_F_A0_data, color='red', label='Ideal Fit')  # This line represents a perfect match

plt.xlabel('Original W/F_A0')
plt.ylabel('Estimated W/F_A0')
plt.title('Estimated vs. Original Values of W/F_A0')
plt.legend()
plt.grid(True)

plt.figure(figsize=(8, 6))
plt.scatter(W_F_A0_data, X_data, color='blue', label='Estimated vs. Original')
plt.show()
