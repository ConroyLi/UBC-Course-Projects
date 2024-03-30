# Define the structure for Gauss-Newton method with Marquardt's modification
import numpy as np
#GNM Method
def objective_function(params, X, Y, model):
    """
    Calculate the objective function (sum of squared residuals).
    
    :param params: numpy array, current parameters
    :param X: numpy array, independent variables
    :param Y: numpy array, dependent variables (experimental values)
    :param model: function, the model predicting Y based on X and params
    :return: float, the value of the objective function
    """
    predictions = model(X, params)
    residuals = Y - predictions
    return np.sum(residuals**2)

def jacobian_matrix(params, X, Y, model, epsilon=1e-6):
    """
    Approximate the Jacobian matrix of the model function with respect to the parameters.
    
    :param params: numpy array, current parameters
    :param X: numpy array, independent variables
    :param Y: numpy array, dependent variables (experimental values)
    :param model: function, the model predicting Y based on X and params
    :param epsilon: float, a small value for numerical differentiation
    :return: numpy array, the Jacobian matrix
    """
    jacobian = np.zeros((len(Y), len(params)))
    for i in range(len(params)):
        params_eps = params.copy()
        params_eps[i] += epsilon
        grad_approx = (model(X, params_eps) - model(X, params)) / epsilon
        jacobian[:, i] = grad_approx
    return jacobian

def gauss_newton_marquardt(X, Y, model, initial_params, max_iterations=1000, tau=1e-3, epsilon=1e-6, convergence_threshold=1e-5):
    """
    Implement the Gauss-Newton method with Marquardt's modification.
    
    :param X: numpy array, independent variables
    :param Y: numpy array, dependent variables (experimental values)
    :param model: function, the model predicting Y based on X and params
    :param initial_params: numpy array, initial guess for the parameters
    :param max_iterations: int, maximum number of iterations
    :param tau: float, damping factor for Marquardt's modification
    :param epsilon: float, a small value for numerical differentiation
    :param convergence_threshold: float, threshold for convergence check
    :return: numpy array, estimated parameters
    """
    params = initial_params.copy()
    for iteration in range(max_iterations):
        J = jacobian_matrix(params, X, Y, model, epsilon)
        residuals = Y - model(X, params)
        H = J.T @ J + tau * np.eye(len(params))  # Hessian approximation with Marquardt's modification
        g = J.T @ residuals
        param_update = np.linalg.inv(H) @ g
        params += param_update
        if np.linalg.norm(param_update) < convergence_threshold:
            print(f"Convergence reached after {iteration+1} iterations.")
            break
    return params

# Question 1

# Model A function definition
def model_A(X, params):
    """
    Predicts the initial rate rAi for Catalytic De-hydrogenation of Sec-butyl Alcohol using Model A.
    
    :param X: numpy array, independent variables (pA - partial pressure of sec-butyl alcohol)
    :param params: numpy array, parameters (kH, kR, KA)
    :return: numpy array, predicted rAi values
    """
    kH, kR, KA = params
    pA = X  # Assuming X is just pA for simplification
    
    R = kH + kH**2/(2*kR) * (1+KA*pA)**2/(KA*pA)
    rAi = R - np.sqrt(R**2 - kH**2)
    return rAi

# Data preparation (simplified example)
pA_values = np.array([1.0, 7.0, 4.0, 10.0, 14.6, 5.5, 8.5, 3.0, 0.22, 1.0])
                     # 1.0, 3.0, 5.0, 7.0, 9.6])  # Pressure (atm)
rAi_values = np.array([0.0392, 0.0416, 0.0416, 0.0326, 0.0247, 0.0415, 0.0376, 0.0420, 0.0295, 0.0410])
                     #  0.0227, 0.0277, 0.0255, 0.0217, 0.0183])  # Initial Rate

# Initial parameters from Table 1 for Model A at 600 oF (for example)
initial_params_A = np.array([7.89e-2, 81.7e-2, 53.5e-2])  # Assuming kH, kR, KA are given in correct units

# Using the gauss_newton_marquardt function to estimate parameters for Model A
estimated_params_A = gauss_newton_marquardt(pA_values, rAi_values, model_A, initial_params_A)

print(estimated_params_A)

# Model B function definition
def model_B(X, params):
    """
    Predicts the initial rate rAi for Catalytic De-hydrogenation of Sec-butyl Alcohol using Model B.
    
    :param X: numpy array, independent variables (pA - partial pressure of sec-butyl alcohol)
    :param params: numpy array, parameters (kH, kR, KA, lambda)
    :return: numpy array, predicted rAi values
    """
    kH, kR, KA= params
    lambda_ = -0.7
    pA = X  # Assuming X is just pA for simplification
    
    rAi = ((kH ** lambda_) + (kR * KA * pA / ((1 + KA * pA)**2))**lambda_) **(1/ lambda_)
    return rAi

# Initial parameters from Table 2 for Model B at 600 oF (as an example)
# Note: Added lambda as -0.7 based on the problem description
initial_params_B = np.array([9.50e-2, 62.8e-2, 51.5e-2])  # Assuming kH, kR, KA are given in correct units

# Using the gauss_newton_marquardt function to estimate parameters for Model B
estimated_params_B = gauss_newton_marquardt(pA_values, rAi_values, model_B, initial_params_B)

print(estimated_params_B)

# Question B
# Define the model function for the Oxidation of Propylene problem
def oxidation_propylene_model(X, params):
    """
    Predicts the rate of propylene disappearance rp.
    
    :param X: numpy array, independent variables [co, cp, n]
    :param params: numpy array, parameters (ko, kp)
    :return: numpy array, predicted rp values
    """
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
initial_params_oxidation = np.array([1300, 0.5])  # Initial guesses for ko and kp

# Use the gauss_newton_marquardt function to estimate parameters for the Oxidation of Propylene
estimated_params_oxidation = gauss_newton_marquardt(X_data, rp_values, oxidation_propylene_model, initial_params_oxidation)

print(estimated_params_oxidation)
