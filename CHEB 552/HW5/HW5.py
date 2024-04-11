# Define the structure for implementing the Gauss-Newton Method with Marquardt's Modification

# Import necessary libraries
# Note: Assuming that numpy, scipy, and matplotlib are available in the coding environment.
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Step 1: Define the ODE system for Toluene Hydrogenation
def toluene_ode_system(concentrations, t, params):
    """
    Defines the system of ODEs for the Toluene Hydrogenation process.
    
    Parameters:
    - concentrations: Array of concentrations [CA, CB, CC]
    - t: Time variable
    - params: Parameters of the system [kH, kD, k2, KA_rel, KC_rel]
    
    Returns:
    - Array of derivatives [dCA/dt, dCB/dt, dCC/dt]
    """
    CA, CB, CC = concentrations
    kH, kD, k2, KA_rel, KC_rel = params
    
    # Rate equations
    r1 = kH * KA_rel * CA / (KA_rel * CA +  CB +  KC_rel*CC)
    r_1 = kD * CB / (KA_rel * CA +  CB +  KC_rel*CC)
    r2 = k2 * CB / (KA_rel * CA +  CB +  KC_rel*CC)
    
    # ODE system
    dCA_dt = -r1 + r_1
    dCB_dt = r1 - r_1 - r2
    dCC_dt = r2
    
    return [dCA_dt, dCB_dt, dCC_dt]

# Example initial parameters and conditions for testing the ODE system
initial_concentrations = [1.0, 0.0, 0.0]  # Initial concentrations of A, B, C
params = [0.023, 0.005, 0.011, 1.9, 1.8]  # Example parameters based on the document
t = [0, 15, 30, 45, 60, 75, 90, 120, 180, 240, 320, 360, 380, 400]  # Time range from 0 to 400 minutes

# Numerical integration of the ODE system
solution = odeint(toluene_ode_system, initial_concentrations, t, args=(params,))
#print(solution)
# Plotting the solution
plt.plot(t, solution[:, 0], label='CA')
plt.plot(t, solution[:, 1], label='CB')
plt.plot(t, solution[:, 2], label='CC')
plt.xlabel('Time (min)')
plt.ylabel('Concentration')
plt.title('Toluene Hydrogenation Kinetics')
plt.legend()
#plt.show()

# Sensitivity Matrix (G) Calculation - Numerical Approximation
def calculate_sensitivity_matrix(system_func, params, concentrations, t, delta=1e-8):
    num_params = len(params)
    num_outputs = len(concentrations)  # Assuming 'concentrations' holds initial values for all outputs
    n = len(t)
    m = num_outputs  # Number of outputs

    # Initialize the sensitivity matrix G with the correct shape: (n*m, p)
    # n * m rows for each time point and output, p columns for each parameter
    sensitivities = np.zeros((n * m, num_params))

    # Calculate the baseline solution with current parameters
    baseline_solution = odeint(system_func, concentrations, t, args=(params,))

    for i in range(num_params):
        # Perturb each parameter by a small delta
        params_perturbed = np.array(params, copy=True)
        params_perturbed[i] += delta
        
        # Calculate the solution with the perturbed parameter
        perturbed_solution, info = odeint(system_func, concentrations, t, args=(params_perturbed,), full_output=1)
        # Calculate the sensitivity for this parameter across all outputs and time points
        # This involves taking the difference between the perturbed and baseline solutions, divided by delta
        # We then reshape this to ensure it's a column in G, representing the sensitivity to parameter i
        sensitivity = (perturbed_solution - baseline_solution) / delta
        sensitivities[:, i] = sensitivity.flatten()  # Flatten to make sure it aligns as a column for each parameter

    return sensitivities


# Iterative Solution Process for the Gauss-Newton Method with Marquardt's Modification
def gauss_newton_marquardt(system_func, initial_params, initial_concentrations, t, experimental_data, max_iterations=100, lambda_factor=10, convergence_threshold=1e-4):
    """
    Implements the Gauss-Newton method with Marquardt's modification for parameter estimation.
    
    Parameters:
    - system_func: The function defining the ODE system.
    - initial_params: Initial guesses for the parameters.
    - initial_concentrations: Initial concentrations for the ODE system.
    - t: Time points for the ODE solution and experimental data.
    - experimental_data: Experimental data for comparison.
    - max_iterations: Maximum number of iterations.
    - lambda_factor: Initial damping factor for Marquardt's modification.
    - convergence_threshold: Convergence criterion based on parameter update magnitude.
    
    Returns:
    - Final parameter estimates.
    """
    params = initial_params
    num_params = len(params)
    for iteration in range(max_iterations):
        # Solve the system with current parameters
        modeled_data = odeint(system_func, initial_concentrations, t, args=(params,))
        
        # Calculate residuals
        residuals = (modeled_data - experimental_data).flatten()
        
        # Calculate sensitivity matrix
        G = calculate_sensitivity_matrix(system_func, params, initial_concentrations, t)
        
        # Compute the Jacobian and the Hessian approximation
        J = np.dot(G.T, residuals)  # Assuming residuals is properly shaped
        H = np.dot(G.T, G) + np.eye(num_params) * lambda_factor

        
        # Parameter update
        param_update = np.linalg.solve(H, J)
        # Assuming param_update should only have as many elements as there are parameters
        params += param_update.flatten()
        
        # Check convergence
        if np.linalg.norm(param_update) < convergence_threshold:
            print(f"Convergence achieved after {iteration + 1} iterations.")
            break
    
    return params
def calculate_confidence_intervals(params, G, residuals):
    """
    Calculates 95% confidence intervals for the estimated parameters.
    
    Parameters:
    - params: Estimated parameters.
    - G: Sensitivity matrix at the optimal parameters.
    - residuals: Differences between the experimental and modeled data.
    
    Returns:
    - Confidence intervals for each parameter.
    """
    
    sigma_squared = np.sum(residuals**2) / (len(residuals) - len(params))
    H_inv = np.linalg.inv(np.dot(G.T, G))  # Approximation of the inverse Hessian
    var_params = sigma_squared * np.diag(H_inv)
    std_errors = np.sqrt(var_params)
    confidence_intervals = [(param - 1.96 * se, param + 1.96 * se) for param, se in zip(params, std_errors)]
    return confidence_intervals

def plot_model_fit(experimental_data, modeled_data, compounds=['CA', 'CB', 'CC']):
    """
    Plots both experimental and modeled data over time for visual comparison.
    
    Parameters:
    - experimental_data: 2D array with time in the first column and experimental concentrations in subsequent columns.
    - modeled_data: 2D array with modeled concentrations for each compound at the times given in experimental_data.
    - compounds: List of compound names corresponding to the columns in experimental_data and modeled_data.
    """
    t = experimental_data[:, 0]  # Extract time points from the first column of experimental_data
    plt.figure(figsize=(12, 8))

    for i, compound in enumerate(compounds, start=1):
        plt.subplot(3, 1, i)
        plt.plot(t, experimental_data[:, i], 'o', label=f'{compound} (Experimental)')
        plt.plot(t, modeled_data[:, i-1], '-', label=f'{compound} (Modeled)')  # Adjust index for modeled_data if necessary
        plt.xlabel('Time (min)')
        plt.ylabel('Concentration')
        plt.title(f'Model Fit for {compound}')
        plt.legend()

    plt.tight_layout()
    plt.show()
experimental_data = np.array([
    [0, 1.000, 0.000, 0.000],
    [15, 0.695, 0.312, 0.001],
    [30, 0.492, 0.430, 0.080],
    [45, 0.276, 0.575, 0.151],
    [60, 0.225, 0.570, 0.195],
    [75, 0.163, 0.575, 0.224],
    [90, 0.134, 0.533, 0.330],
    [120, 0.064, 0.462, 0.471],
    [180, 0.056, 0.362, 0.580],
    [240, 0.041, 0.211, 0.747],
    [320, 0.031, 0.146, 0.822],
    [360, 0.022, 0.080, 0.898],
    [380, 0.021, 0.070, 0.909],
    [400, 0.019, 0.073, 0.908]
])
# Extract time points and concentrations from the experimental data for analysis
time_points = experimental_data[:, 0]
experimental_concentrations = experimental_data[:, 1:]

# Initial parameter guesses for the optimization process
initial_params = [0.023, 0.005, 0.011, 1.9, 1.8]  # As reported in the problem statement

# Perform the optimization to estimate the parameters
estimated_params = gauss_newton_marquardt(
    system_func=toluene_ode_system,
    initial_params=initial_params,
    initial_concentrations=experimental_concentrations[0],
    t=time_points,
    experimental_data=experimental_concentrations
)

# Use the estimated parameters to generate model predictions
modeled_data = odeint(toluene_ode_system, experimental_concentrations[0], time_points, args=(estimated_params,))

# Visualize the fit of the model to the experimental data
plot_model_fit(experimental_data, modeled_data)

# Calculate and print the confidence intervals for the estimated parameters
confidence_intervals = calculate_confidence_intervals(estimated_params, G, residuals)
print("Confidence Intervals for Parameters:", confidence_intervals)

