import numpy as np
from numpy import log
from scipy.stats import t

# Assuming these are your original data points
p_data = np.array([0.3427, 0.5134, 0.6837, 0.8535, 1.0229, 1.1918, 1.5284, 1.7821, 2.1239, 2.5752, 3.2134, 3.5572, 3.679, 4.39, 4.8735, 5.1615])
q_data = np.array([1.9267, 2.1711, 2.3449, 2.453, 2.5799, 2.6612, 2.7937, 2.851, 2.9159, 2.9991, 3.1062, 3.1438, 3.1551, 3.273, 3.3012, 3.3228])

# Model function with parameter scaling
def model(p, params_scaled, scale_factors):
    params = params_scaled * scale_factors
    q_sat, k, t = params
    return (q_sat * k * p) / ((1 + (k * p)**t)**(1/t))

# Residual function for scaled parameters
def residual(params_scaled, p, q_obs, scale_factors):
    return q_obs - model(p, params_scaled, scale_factors)

# Jacobian function adjusted for scaled parameters (simplified for demonstration)
def jacobian(p, params_scaled, scale_factors):
    q_sat, k, t = params_scaled
    # Implement the Jacobian calculation for the scaled parameters here
    # This is a placeholder; you'll need to adjust this based on your actual Jacobian
    J = np.zeros((len(p), len(params_scaled)))
    # Remember to adjust the computation for the scaling
    J[:, 0] = k * p / (1 + (k * p) ** t)**(1/t)

    # Partial derivative with respect to k
    J[:, 1] = p*q_sat*(k**t*p**t + 1)**(-1 - 1/t)

    # Partial derivative with respect to t
    J[:, 2] =  -log(((k*p)**(t*(k*p)**t)*((k*p)**t + 1)**(-(k*p)**t - 1))**(k*p*q_sat/(t**2*((k*p)**t + 1)**((t + 1)/t))))

    return J

# Levenberg-Marquardt algorithm adapted for parameter scaling
def levenberg_marquardt(p, q, params_guess, scale_factors, max_iter=100):
    params_scaled = params_guess / scale_factors
    
    for i in range(max_iter):
        # Here, implement the optimization steps using scaled parameters
        # This is a placeholder loop. You'll need to fill in the optimization details
        pass
    
    # After optimization, scale the parameters back
    estimated_params = params_scaled * scale_factors
    
    # Covariance matrix calculation (adjusted for scaled parameters)
    # Placeholder for covariance matrix calculation
    covar = np.eye(len(params_scaled))  # Placeholder, replace with actual calculation
    
    # Calculate standard errors and confidence intervals
    se_scaled = np.sqrt(np.diag(covar))
    se = se_scaled / scale_factors  # Adjust SEs based on scaling
    
    df = len(p) - len(params_scaled)
    t_critical = t.ppf(0.975, df)
    ci_lower = estimated_params - t_critical * se
    ci_upper = estimated_params + t_critical * se
    
    return estimated_params, np.vstack((ci_lower, ci_upper))

# Initial guesses and scaling factors
initial_params = np.array([5, 18, 0.7])
scale_factors = np.array([1, 10, 0.1])  # Example scaling, adjust as needed

# Perform optimization with scaling
estimated_params, confidence_intervals = levenberg_marquardt(p_data, q_data, initial_params, scale_factors)

print("Estimated Parameters:", estimated_params)
print("Confidence Intervals:\n", confidence_intervals)
