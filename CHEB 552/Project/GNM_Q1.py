#0.3427, 0.5134, 0.6837, 0.8535, 1.0229, 1.1918, 1.5284, 1.7821, 2.1239, 2.5752, 3.2134, 3.5572, 3.679, 4.39, 4.8735, 5.1615
#1.9267, 2.1711, 2.3449, 2.453, 2.5799, 2.6612, 2.7937, 2.851, 2.9159, 2.9991, 3.1062, 3.1438, 3.1551, 3.273, 3.3012, 3.3228
# Let's start by implementing the model, residual and Jacobian functions and the Levenberg-Marquardt algorithm

# Let's carefully check the Jacobian implementation and the parameter update step
# and make the necessary corrections.
from numpy import log
import numpy as np
from scipy.stats import t
import matplotlib.pyplot as plt
def model(p, q_sat, k, t):
    # Ensure all operations are element-wise
    p = np.array(p, dtype=float)
    term = (1 + (k * p) ** t)
    return q_sat * k * p / term ** (1 / t)
def residual(params, p, q_obs):
    q_sat, k, t = params
    return q_obs - model(p, q_sat, k, t)
def jacobian(p, params):
    q_sat, k, t = params
    # Ensure p is a numpy array for element-wise operations
    p = np.array(p, dtype=float)
    J = np.zeros((len(p), len(params)))

    # Compute the model to avoid recomputing
    # q_pred = model(p, q_sat, k, t)

    # Partial derivative with respect to q_sat
    J[:, 0] = k * p / (1 + (k * p) ** t)**(1/t)

    # Partial derivative with respect to k
    J[:, 1] = p*q_sat*(k**t*p**t + 1)**(-1 - 1/t)

    # Partial derivative with respect to t
    J[:, 2] =  -log(((k*p)**(t*(k*p)**t)*((k*p)**t + 1)**(-(k*p)**t - 1))**(k*p*q_sat/(t**2*((k*p)**t + 1)**((t + 1)/t))))

    return J

# Adjust levenberg_marquardt function to correct the update rule and lambda adjustments

def levenberg_marquardt(p, q, params_guess, max_iter=100, lambda_inc=10, lambda_dec=10, lambda_init=0.01):
    params = np.array(params_guess, dtype=float)
    n = len(q)
    iteration_details = []

    # Initial calculation
    r = residual(params, p, q)
    J = jacobian(p, params)
    
    A = J.T @ J
    # print(A)
    g = J.T @ r
    lambda_ = lambda_init * np.max(np.diag(A))
    objective = np.sum(r**2)

    for i in range(max_iter):
        # Solve for parameter update
        try:
            delta = np.linalg.solve(A + lambda_ * np.eye(len(params)), g)
        except np.linalg.LinAlgError:
            print(f"Iteration {i}: Singular matrix encountered. Adjusting lambda.")
            lambda_ *= lambda_inc
            continue

        # Check if improvement
        new_params = params + delta
        new_r = residual(new_params, p, q)
        new_objective = np.sum(new_r**2)

        if new_objective < objective:
            # Update is successful, decrease lambda (move towards Gauss-Newton method)
            lambda_ /= lambda_dec
            params = new_params
            objective = new_objective
            r = new_r
            J = jacobian(p, params)
            A = J.T @ J
            g = J.T @ r
            #print(f"Iteration {i}: Success - Parameters updated to: {params} with objective: {objective}")
        else:
            # Update is not successful, increase lambda (move towards gradient descent)
            lambda_ *= lambda_inc
            #print(f"Iteration {i}: No improvement - Lambda adjusted to: {lambda_}")

        # Store iteration details
        iteration_details.append((i, objective, params))

        # Termination condition on gradient (close to zero gradient)
        if np.linalg.norm(g, np.inf) < 1e-4:
            # print("Termination condition met: Gradient close to zero.")
            break

    # Calculate approximate covariance matrix
    covar = np.linalg.pinv(J.T @ J) 
    # print(covar)
    # Compute the 95% confidence intervals, assuming large sample size for degrees of freedom
    se = np.sqrt(np.diag(covar))
    n = len(p)
    pp = len(params_guess)
    df = n - pp

# Calculate the t-score for 95% confidence interval
    t_score = t.ppf(0.975, df)

# Replace the z-score with the t-score in the confidence interval calculation
    ci_95 = t_score * se/100

    return params, iteration_details, ci_95
# Initial guesses for the parameters
initial_params_303 = [4, 17, 0.6]
initial_params_338 = [6, 9, 0.5]
initial_params_373 = [5, 2, 0.4]
initial_params_423 = [5.5, 0.7, 0.3]
initial_params_473 = [3, 0.3, 0.3]

p_data_303 = [0, 0.3427, 0.5134, 0.6837, 0.8535, 1.0229, 1.1918, 1.5284, 1.7821, 2.1239, 2.5752, 3.2134, 3.5572, 3.679, 4.39, 4.8735, 5.1615]
q_data_303 = [0, 1.9267, 2.1711, 2.3449, 2.453, 2.5799, 2.6612, 2.7937, 2.851, 2.9159, 2.9991, 3.1062, 3.1438, 3.1551, 3.273, 3.3012, 3.3228]

p_data_338 = [0, 0.3559, 0.5332, 0.71, 0.8863, 1.0622, 1.2377, 1.4126, 1.5872, 1.8313, 2.1936, 2.7723, 3.2543, 3.6185, 3.759, 4.5065, 5.004, 5.2772]
q_data_338 = [0, 1.2594, 1.4662, 1.6353, 1.7556, 1.8581, 1.9577, 2.0324, 2.1109, 2.2044, 2.3167, 2.4243, 2.4803, 2.5766, 2.6024, 2.7204, 2.7979, 2.813]

p_data_373 = [0, 0.3651, 0.5470, 0.7284, 0.9093, 1.0898, 1.2697, 1.4493, 1.6283, 1.8663, 2.2458, 2.8411, 3.2951, 3.6798, 3.8230, 4.6036, 5.0910, 5.3930]
q_data_373 = [0, 0.5874, 0.7237, 0.8365, 0.9305, 0.9803, 1.0818, 1.1316, 1.2030, 1.2768, 1.3440, 1.4685, 1.5268, 1.6034, 1.6198, 1.7199, 1.7857, 1.8120]

p_data_423 = [0, 0.3849, 0.5766, 0.7678, 0.9586, 1.1488, 1.3385, 1.5277, 1.7165, 1.9576, 2.3329, 2.9556, 3.4313, 3.8025, 3.9350, 4.7008, 5.2216, 5.5087]
q_data_423 = [0, 0.3102, 0.3853, 0.4525, 0.4944, 0.5503, 0.5954, 0.6297, 0.6767, 0.7265, 0.7810, 0.8623, 0.9253, 0.9662, 0.9826, 1.0526, 1.0982, 1.1194]

p_data_473 = [0, 0.3981, 0.5964, 0.7941, 0.9914, 1.1881, 1.3843, 1.5801, 1.7753, 2.0347, 2.4199, 3.0243, 3.5130, 3.8945, 4.0469, 4.7785, 5.2651, 5.5550]
q_data_473 = [0, 0.1175, 0.1551, 0.1833, 0.2054, 0.2162, 0.2547, 0.2693, 0.2923, 0.3163, 0.3407, 0.3839, 0.4145, 0.4370, 0.4474, 0.4737, 0.5066, 0.5202]

# Uncomment the function call to execute the algorithm with the data
estimated_params_303, iterations_info_303, confidence_intervals_303 = levenberg_marquardt(p_data_303, q_data_303, initial_params_303)
estimated_params_338, iterations_info_338, confidence_intervals_338 = levenberg_marquardt(p_data_338, q_data_338, initial_params_338)
estimated_params_373, iterations_info_373, confidence_intervals_373 = levenberg_marquardt(p_data_373, q_data_373, initial_params_373)
estimated_params_423, iterations_info_423, confidence_intervals_423 = levenberg_marquardt(p_data_423, q_data_423, initial_params_423)
estimated_params_473, iterations_info_473, confidence_intervals_473 = levenberg_marquardt(p_data_473, q_data_473, initial_params_473)
'''
print(np.round(estimated_params_303, 2))
print(np.round(confidence_intervals_303, 2))

print(np.round(estimated_params_338, 2))
print(np.round(confidence_intervals_338, 2))

print(np.round(estimated_params_373, 2))
print(np.round(confidence_intervals_373, 2))

print(np.round(estimated_params_423, 2))
print(np.round(confidence_intervals_423, 2))

print(np.round(estimated_params_473, 2))
print(np.round(confidence_intervals_473, 2))
'''
''' For GNM
q_est_303 = model(p_data_303,estimated_params_303[0],estimated_params_303[1],estimated_params_303[2])
q_est_338 = model(p_data_338,estimated_params_338[0],estimated_params_338[1],estimated_params_338[2])
q_est_373 = model(p_data_373,estimated_params_373[0],estimated_params_373[1],estimated_params_373[2])
q_est_423 = model(p_data_423,estimated_params_423[0],estimated_params_423[1],estimated_params_423[2])
q_est_473 = model(p_data_473,estimated_params_473[0],estimated_params_473[1],estimated_params_473[2])
#print(q_est_303)
#print(q_est_338)

'''
'''For NMA'''
q_est_303 = model(p_data_303,4.31,17.81,0.45)
q_est_338 = model(p_data_338,5.58,9.99,0.34)
q_est_373 = model(p_data_373,5.68,2.1,0.33)
q_est_423 = model(p_data_423,5.77,0.74,0.31)
q_est_473 = model(p_data_473,2.92,0.3,0.35)
#print(q_est_303)
#print(q_est_338)


plt.figure(figsize=(10, 18))

#plt.title('Experiment Data vs. Estimated Data of q for Problem 1')
plt.subplot(5,1,1)
plt.title('Experiment Data vs. Estimated Data of q for Problem 1 T=303K')
plt.scatter(p_data_303, q_data_303, color='blue', label='Experiment Data')
plt.plot(p_data_303,q_est_303, color='red', label='Estimated Data') 
plt.legend()
plt.grid(True) 
plt.ylabel('q')
plt.subplot(5,1,2)
plt.title('T=338K')
plt.scatter(p_data_338, q_data_338, color='blue', label='Experiment Data')
plt.plot(p_data_338,q_est_338, color='red', label='Estimated Data') 
plt.legend()
plt.grid(True)
plt.ylabel('q')
plt.subplot(5,1,3)
plt.title('T=373K')
plt.scatter(p_data_373, q_data_373, color='blue', label='Experiment Data')
plt.plot(p_data_373,q_est_373, color='red', label='Estimated Data') 
plt.legend()
plt.grid(True)
plt.ylabel('q')


plt.subplot(5,1,4)
plt.title('T=423K')
plt.scatter(p_data_423, q_data_423, color='blue', label='Experiment Data')
plt.plot(p_data_423,q_est_423, color='red', label='Estimated Data') 
plt.legend()
plt.grid(True)
plt.ylabel('q')
plt.subplot(5,1,5)
plt.title('T=473K')
plt.scatter(p_data_473, q_data_473, color='blue', label='Experiment Data')
plt.plot(p_data_473,q_est_473, color='red', label='Estimated Data') 
plt.legend()
plt.grid(True)
plt.xlabel('p')
plt.ylabel('q')

plt.show()

