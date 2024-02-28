import numpy as np
from scipy.optimize import minimize

# Given values
P = 750
w = np.array([-10.021, -21.096, -37.986, -9.846, -28.653, -18.918, -28.032, -14.640, -30.594, -26.111])

# Objective function
def objective(x):
    total_x = np.sum(x)
    return np.sum(x * (w + np.log(P) + np.log(x / total_x)))

# Constraints
constraints = [
    {'type': 'eq', 'fun': lambda x: x[0] + 2*x[1] + 2*x[2] + x[5] + x[9] - 2},
    {'type': 'eq', 'fun': lambda x: x[3] + 2*x[4] + x[5] + x[6] - 1},
    {'type': 'eq', 'fun': lambda x: x[2] + x[6] + x[7] + 2*x[8] + x[9] - 1}
]

# Initial guess
x0 = np.ones(10) / 10  # Example initial guess

# Solve the optimization problem
result = minimize(objective, x0, method='SLSQP', constraints=constraints)

print("Optimal x:", result.x)
print("Minimum f(x):", result.fun)
