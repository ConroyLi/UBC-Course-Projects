from sympy import symbols, log, solve, Eq
import numpy as np
from scipy.optimize import minimize

# Given values
P = 750
w = [-10.021, -21.096, -37.986, -9.846, -28.653, -18.918, -28.032, -14.640, -30.594, -26.111]

# Define symbols for xi
x = symbols('x1:11')  # x1 through x10

# Objective function (for symbolic setup, actual optimization will use a numerical approach)
f_x = sum(x[i] * (w[i] + log(P) + log(x[i]/sum(x))) for i in range(10))

# Material balance constraints
constraints = [
    Eq(x[0] + 2*x[1] + 2*x[2] + x[5] + x[9], 2),
    Eq(x[3] + 2*x[4] + x[5] + x[6], 1),
    Eq(x[2] + x[6] + x[7] + 2*x[8] + x[9], 1)
]

# Due to the complexity of directly solving this symbolic representation, 
# we will switch to a numerical optimization method using scipy.optimize.

# Define the numerical objective function
def numerical_obj(x):
    total_x = np.sum(x)
    return np.sum(x * (np.array(w) + np.log(P) + np.log(x/total_x)))

# Initial guess (based on constraints, a simple guess could be an equal distribution)
x0 = np.full(10, 0.1)  # Starting point for optimization

# Define the constraints in a form suitable for scipy.optimize
cons = [
    {'type': 'eq', 'fun': lambda x: x[0] + 2*x[1] + 2*x[2] + x[5] + x[9] - 2},
    {'type': 'eq', 'fun': lambda x: x[3] + 2*x[4] + x[5] + x[6] - 1},
    {'type': 'eq', 'fun': lambda x: x[2] + x[6] + x[7] + 2*x[8] + x[9] - 1}
]

# Bounds to ensure non-negative solutions
bounds = [(0, None) for _ in range(10)]

# Perform the optimization
result = minimize(numerical_obj, x0, method='SLSQP', constraints=cons, bounds=bounds)

result.x, result.fun
