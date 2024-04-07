from sympy import symbols, diff, simplify

# Define the symbols
q_sat, k, p, t = symbols('q_sat k p t', real=True, positive=True)

# Define the model equation
q_model = (q_sat * k * p) / ((1 + (k * p)**t)**(1/t))

# Compute the partial derivatives symbolically
partial_q_sat = simplify(diff(q_model, q_sat))
partial_k = simplify(diff(q_model, k))
partial_t = simplify(diff(q_model, t))

# Print the results
print("Partial derivative with respect to q_sat:", partial_q_sat)
print("Partial derivative with respect to k:", partial_k)
print("Partial derivative with respect to t:", partial_t)
