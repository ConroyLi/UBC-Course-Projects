import numpy as np
import matplotlib.pyplot as plt

# Constants
mu = 1  # Shear modulus, arbitrary units
b = 1   # Burgers vector magnitude, arbitrary units

R = np.linspace(0.1, 5, 500)  # Radii from 0.1 to 5 units

# Calculate the strain energy per unit length U as a function of R
U = (mu * b**2 * R**2) / (16 * np.pi)

# Plotting the strain energy per unit length
plt.figure(figsize=(8, 5))
plt.plot(R, U, label='$U = \\frac{\\mu b^2 R^2}{16 \\pi}$')
plt.xlabel('Radius $R$')
plt.ylabel('Strain Energy per Unit Length $U$')
plt.title('Strain Energy per Unit Length vs Radius')
plt.legend()
plt.grid(True)
plt.show()
