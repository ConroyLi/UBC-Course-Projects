# Define the structure for implementing the Gauss-Newton Method with Marquardt's Modification

# Import necessary libraries
# Note: Assuming that numpy, scipy, and matplotlib are available in the coding environment.
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Step 1: Define the ODE system for Toluene Hydrogenation
def methylester_ode(t, y, params):
    CA, CB, CC, CD = y
    k1, k2, k3, KB_rel, KC_rel, KD_rel = params
    dom = CA + KB_rel*CB + KC_rel*CC + KD_rel*CD
    r1 = k1*CA/dom
    r2 = k2*KB_rel*CB/dom
    r3 = k3*KC_rel*CC/dom

    dCA_dt = -r1
    dCB_dt = r1-r2
    dCC_dt = r2-r3
    dCD_dt = r3

    return [dCA_dt, dCB_dt, dCC_dt, dCD_dt]


# Initial parameters and conditions for testing the ODE system
initial_concentrations = [0.1012, 0.221, 0.657, 0.0208]  # Initial concentrations of A, B, C
params = [3, 0.1, 0.001, 50, 0.1, 5]  # Example parameters based on the document
t = [0,10,14,19,24,34,69,124]  # Time range from 0 to 400 minutes
t_span = [t[0], t[-1]]
# Numerical integration of the ODE system
sol = solve_ivp(methylester_ode, t_span, initial_concentrations, args=(params,), t_eval=t, method='BDF')
solution = sol.y.T
#print(solution)
# Plotting the solution
plt.plot(t, solution[:, 0], label='CA')
plt.plot(t, solution[:, 1], label='CB')
plt.plot(t, solution[:, 2], label='CC')
plt.plot(t, solution[:, 3], label='CD')
plt.xlabel('Time (min)')
plt.ylabel('Concentration')
plt.title('Hydrogenation of Methylesters')
plt.legend()
#plt.show()

# Sensitivity Matrix (G) Calculation - Numerical Approximation
def calculate_sensitivity_matrix(system_func, params, concentrations, t, delta=1e-6):
    num_params = len(params)
    num_outputs = len(concentrations)  # Assuming 'concentrations' holds initial values for all outputs
    n = len(t)
    m = num_outputs  # Number of outputs

    # Initialize the sensitivity matrix G with the correct shape: (n*m, p)
    # n * m rows for each time point and output, p columns for each parameter
    sensitivities = np.zeros((n * m, num_params))

    # Calculate the baseline solution with current parameters
    sol = solve_ivp(system_func, t_span, initial_concentrations, args=(params,), t_eval=t, method='BDF')
    solution = sol.y.T

    for i in range(num_params):
        # Perturb each parameter by a small delta
        params_perturbed = np.array(params, copy=True)
        params_perturbed[i] += delta
        
        # Calculate the solution with the perturbed parameter
        perturbed_solution = solve_ivp(system_func, t_span, initial_concentrations, args=(params_perturbed,), t_eval=t, method='BDF')
        # Calculate the sensitivity for this parameter across all outputs and time points
        # This involves taking the difference between the perturbed and baseline solutions, divided by delta
        # We then reshape this to ensure it's a column in G, representing the sensitivity to parameter i
        sensitivity = (perturbed_solution.y.T - solution) / delta
        sensitivities[:, i] = sensitivity.flatten()  # Flatten to make sure it aligns as a column for each parameter

    return sensitivities
# Iterative Solution Process for the Gauss-Newton Method with Marquardt's Modification
def gauss_newton_marquardt(system_func, initial_params, initial_concentrations, t, experimental_data, max_iterations=10, lambda_factor=1, convergence_threshold=1e-4):
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
    iteration_details = []
    for iteration in range(max_iterations):
        # Solve the system with current parameters
        modeled_data = solve_ivp(system_func, t_span, initial_concentrations, args=(params,), t_eval=t, method='BDF')
        modeled_data = modeled_data.y.T
        #print(modeled_data)
        # Calculate residuals
        residuals = (modeled_data - experimental_data).flatten()
        objective_value = np.sum(residuals**2)
        # Calculate sensitivity matrix
        G = calculate_sensitivity_matrix(system_func, params, initial_concentrations, t)
        
        # Compute the Jacobian and the Hessian approximation
        J = np.dot(G.T, residuals)  # Assuming residuals is properly shaped
        H = np.dot(G.T, G) + np.eye(num_params) * lambda_factor

        
        # Parameter update
        param_update = np.linalg.solve(H, J)
        for i in range(len(params)):
            if params[i] + param_update[i] < 0:
                param_update[i] = 0.1* params[i]
        # Assuming param_update should only have as many elements as there are parameters
        params += param_update.flatten()
        #print(params)
        #print(residuals)
        #params = np.maximum(params, 1e-6)
        # Check convergence
        iteration_details.append({
            "iteration": iteration,
            "objective_value": objective_value,
            "parameters": params.copy()  # Store a copy of the parameter array
        })
        if np.linalg.norm(param_update) < convergence_threshold:
            print(f"Convergence achieved after {iteration + 1} iterations.")
            break
    sigma_squared = np.sum(residuals**2) / (len(residuals) - len(params))
    H_inv = np.linalg.pinv(np.dot(G.T, G))  # Approximation of the inverse Hessian
    var_params = sigma_squared * np.diag(H_inv)
    std_errors = np.sqrt(var_params)
    confidence_intervals = [(param - 1.96 * se, param + 1.96 * se) for param, se in zip(params, std_errors)]
    return params, iteration_details, confidence_intervals

    
    
   

def plot_model_fit(experimental_data, modeled_data, compounds=['CA', 'CB', 'CC',"CD"]):
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
        plt.subplot(4, 1, i)
        plt.plot(t, experimental_data[:, i], 'o', label=f'{compound} (Experimental)')
        plt.plot(t, modeled_data[:, i-1], '-', label=f'{compound} (Modeled)')  # Adjust index for modeled_data if necessary
        plt.xlabel('Time (min)')
        plt.ylabel('Concentration')
        plt.title(f'Model Fit for {compound}')
        plt.legend()

    plt.tight_layout()
    plt.show()
experimental_data = np.array([
    [0, 0.1012, 0.2210, 0.6570, 0.0208],
    [10, 0.0150, 0.1064, 0.6941, 0.1977],
    [14, 0.0044, 0.0488, 0.6386, 0.3058],
    [19, 0.0028, 0.0242, 0.5361, 0.4444],
    [24, 0.0029, 0.0015, 0.3956, 0.6055],
    [34, 0.0017, 0.0005, 0.2188, 0.7808],
    [69, 0.0003, 0.0004, 0.0299, 0.9680],
    [124, 0.0001, 0.0002, 0.0001, 0.9982]
])
# Extract time points and concentrations from the experimental data for analysis
time_points = experimental_data[:, 0]
experimental_concentrations = experimental_data[:, 1:]

# Initial parameter guesses for the optimization process
initial_params = [3, 0.1, 0.001, 50, 0.1, 5]  # As reported in the problem statement
#[0.023, 0.005, 0.011, 1.9, 1.8]
# Perform the optimization to estimate the parameters
estimated_params, iteration_details, CI = gauss_newton_marquardt(
    system_func=methylester_ode,
    initial_params=initial_params,
    initial_concentrations=experimental_concentrations[0],
    t=time_points,
    experimental_data=experimental_concentrations
)

# Use the estimated parameters to generate model predictions
modeled_data = solve_ivp(methylester_ode, t_span, initial_concentrations, args=(estimated_params,), t_eval=t, method='BDF')
modeled_data = modeled_data.y.T
#print(estimated_params)
#print(modeled_data)
# Visualize the fit of the model to the experimental data
plot_model_fit(experimental_data, solution)

# Calculate and print the confidence intervals for the estimated parameters


def print_latex_table(iteration_details):
    # Start the table, define the alignment for each column
    latex_table = "\\begin{table}[ht]\n\\centering\n"
    latex_table += "\\begin{tabular}{|c|c|c|}\n\\hline\n"
    
    # Adding the header row
    latex_table += "Iteration & Objective Value & Parameters \\\\ \\hline\n"
    
    # Iterate through the details, adding each set of details to the table
    for detail in iteration_details:
        iteration = detail["iteration"]
        obj_value = detail["objective_value"]
        # Convert parameters to list if it's a NumPy array
        if isinstance(detail["parameters"], np.ndarray):
            params = detail["parameters"].tolist()
        else:
            params = detail["parameters"]
        params_formatted = ', '.join(f"{p:.4f}" for p in params)  # Format parameters with 4 decimal places
        
        # Add row to table
        latex_table += f"{iteration} & {obj_value:.4f} & {params_formatted} \\\\ \\hline\n"
    
    # End the table
    latex_table += "\\end{tabular}\n"
    latex_table += "\\caption{Iteration Details of the Gauss-Newton Marquardt Optimization}\n"
    latex_table += "\\label{tab:iteration_details}\n"
    latex_table += "\\end{table}"
    
    print(latex_table)

# Example usage, assuming `iteration_details` is the list of dictionaries collected from your optimization function
print_latex_table(iteration_details)
#print(iteration_details)
#CI = ", ".join("(%0.3f, %0.3f)" % (low, high) for low, high in CI)
#print("Confidence Intervals for Parameters:", CI)