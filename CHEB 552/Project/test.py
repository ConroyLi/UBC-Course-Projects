from sympy import symbols, diff, simplify

# Define the symbols
q_sat, k, p, t = symbols('q_sat k p t', real=True, positive=True)

# Define the model equation
q_model = (q_sat * k * p) / ((1 + (k * p)**t)**(1/t))

# Compute the partial derivatives symbolically
partial_q_sat = simplify(diff(q_model, q_sat))
partial_k = simplify(diff(q_model, k))
partial_t = simplify(diff(q_model, t))
'''
# Print the results
print("Partial derivative with respect to q_sat:", partial_q_sat)
print("Partial derivative with respect to k:", partial_k)
print("Partial derivative with respect to t:", partial_t)
'''
import sympy as sp

# Define the symbolic variables
X = sp.symbols('X')
k1, k2, k3 = sp.symbols('k1 k2 k3')

# Define the constants as symbols
HOAC = sp.symbols('[HOAC]')
O2 = sp.symbols('[O2]')

# Define the model function symbolically
model_func =  ((1 + k3 * O2)**2 * (k2 * HOAC* X - sp.log(1-X)) / (k1 * HOAC * O2))

# Compute the symbolic derivatives of the model function with respect to each parameter
derivative_k1 = model_func.diff(k1)
derivative_k2 = model_func.diff(k2)
derivative_k3 = model_func.diff(k3)

# Convert the symbolic derivatives to numerical functions for evaluation
func_derivative_k1 = sp.lambdify((X, k1, k2, k3, HOAC, O2), derivative_k1, 'numpy')
func_derivative_k2 = sp.lambdify((X, k1, k2, k3, HOAC, O2), derivative_k2, 'numpy')
func_derivative_k3 = sp.lambdify((X, k1, k2, k3, HOAC, O2), derivative_k3, 'numpy')

print(derivative_k1)
print(derivative_k2)
print(derivative_k3)
