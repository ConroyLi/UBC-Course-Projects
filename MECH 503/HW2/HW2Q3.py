import numpy as np
import matplotlib.pyplot as plt

# Given values in the problem
theta = 0.1  # Radians should be degree
lambda_ = 0.001
f = 1.001

# Define the matrices F1, F2, F3, F4, F5
F1 = np.array([[1 - lambda_, 0, 0],
               [0, 1 - lambda_, 0],
               [0, 0, 1 - lambda_]])

F2 = np.array([[1 + lambda_, 0, 0],
               [0, 1 + lambda_, 0],
               [0, 0, 1 + lambda_]])

F3 = np.array([[np.cos(theta), -np.sin(theta), 0],
               [np.sin(theta), np.cos(theta), 0],
               [0, 0, 1]])

F4 = np.array([[1.001, 0.002, -0.003],
               [0.003, 1.002, -0.001],
               [0.002, -0.002, 1.]])

F5 = np.array([[f * np.cos(theta), -f * np.sin(theta), 0],
               [f * np.sin(theta), f * np.cos(theta), 0],
               [0, 0, 1 + lambda_]])

# Compute the Jacobian of each F
J_F1 = np.linalg.det(F1)
J_F2 = np.linalg.det(F2)
J_F3 = np.linalg.det(F3)
J_F4 = np.linalg.det(F4)
J_F5 = np.linalg.det(F5)
print(J_F1,J_F2,J_F3,J_F4,J_F5)

# Define points that form a square and a line
p0 = np.array([0, 0, 0])
p1 = np.array([0, 1, 0])
p2 = np.array([1, 1, 0])
p3 = np.array([1, 0, 0])
line_points = np.array([p0, p1]).T
square_points = np.array([p0, p1, p2, p3, p0]).T

# Function to apply transformation and plot
def plot_transformation(F, points, transformation_name):
    transformed_points = F @ points
    plt.plot(transformed_points[0, :], transformed_points[1, :], label=transformation_name)
    plt.legend()
    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)
    return transformed_points

'''
# Apply transformations to the line and square
plt.subplot(5, 2, 1)
transformed_line_F1 = plot_transformation(F1, line_points, 'F1 on line')


plt.subplot(5, 2, 3)
transformed_line_F2 = plot_transformation(F2, line_points, 'F2 on line')

plt.subplot(5, 2, 5)
transformed_line_F3 = plot_transformation(F3, line_points, 'F3 on line')

plt.subplot(5, 2, 7)
transformed_line_F4 = plot_transformation(F4, line_points, 'F4 on line')

plt.subplot(5, 2, 9)
transformed_line_F5 = plot_transformation(F5, line_points, 'F5 on line')

plt.subplot(5, 2, 2)
transformed_line_F1 = plot_transformation(F1, square_points, 'F1 on  square')

plt.subplot(5, 2, 4)
transformed_line_F2 = plot_transformation(F2, square_points, 'F2 on  square')

plt.subplot(5, 2, 6)
transformed_line_F3 = plot_transformation(F3, square_points, 'F3 on  square')

plt.subplot(5, 2, 8)
transformed_line_F4 = plot_transformation(F4, square_points, 'F4 on  square')

plt.subplot(5, 2, 10)
transformed_line_F5 = plot_transformation(F5, square_points, 'F5 on  square')

'''
'''
F23 = F2@F3
plt.subplot(2, 1, 1)
transformed_line_F23 = plot_transformation(F23, square_points, 'F23 on  square')
# F12 = F1@F2
plt.subplot(2, 1, 2)
transformed_line_F5 = plot_transformation(F5, square_points, 'F5 on  square')

#print(F1,F2,F3,F23)
'''


plt.show()

# For Q4 define two fiber elements
fiber_f1 = np.array([[0, 0], [0, 1], [0, 0]])  # from p0 to p1
fiber_f2 = np.array([[0, 1], [0, 0], [0, 0]])  # from p0 to p3

# Calculate the change in angle between fibers f1 and f2
# Calculate the change in length for fibers f1 and f2 using F4 and F5
transformed_fiber_f1_F4 = F4 @ fiber_f1
transformed_fiber_f2_F4 = F4 @ fiber_f2
transformed_fiber_f1_F5 = F5 @ fiber_f1
transformed_fiber_f2_F5 = F5 @ fiber_f2

# Original and transformed lengths of f1 and f2
original_length_f1 = np.linalg.norm(fiber_f1[:, 1] - fiber_f1[:, 0])
original_length_f2 = np.linalg.norm(fiber_f2[:, 1] - fiber_f2[:, 0])
length_f1_F4 = np.linalg.norm(transformed_fiber_f1_F4[:, 1] - transformed_fiber_f1_F4[:, 0])
length_f2_F4 = np.linalg.norm(transformed_fiber_f2_F4[:, 1] - transformed_fiber_f2_F4[:, 0])
length_f1_F5 = np.linalg.norm(transformed_fiber_f1_F5[:, 1] - transformed_fiber_f1_F5[:, 0])
length_f2_F5 = np.linalg.norm(transformed_fiber_f2_F5[:, 1] - transformed_fiber_f2_F5[:, 0])

# For angle calculation, we consider the vectors formed by f1 and f2
vector_f1 = fiber_f1[:, 1] - fiber_f1[:, 0]
vector_f2 = fiber_f2[:, 1] - fiber_f2[:, 0]

transformed_vector_f1_F4 = transformed_fiber_f1_F4[:, 1] - transformed_fiber_f1_F4[:, 0]
transformed_vector_f2_F4 = transformed_fiber_f2_F4[:, 1] - transformed_fiber_f2_F4[:, 0]
transformed_vector_f1_F5 = transformed_fiber_f1_F5[:, 1] - transformed_fiber_f1_F5[:, 0]
transformed_vector_f2_F5 = transformed_fiber_f2_F5[:, 1] - transformed_fiber_f2_F5[:, 0]

# Compute angles using dot product
angle_f1_f2_original = np.arccos(np.dot(vector_f1, vector_f2) / (np.linalg.norm(vector_f1) * np.linalg.norm(vector_f2)))
angle_f1_f2_F4 = np.arccos(np.dot(transformed_vector_f1_F4, transformed_vector_f2_F4) / (np.linalg.norm(transformed_vector_f1_F4) * np.linalg.norm(transformed_vector_f2_F4)))
angle_f1_f2_F5 = np.arccos(np.dot(transformed_vector_f1_F5, transformed_vector_f2_F5) / (np.linalg.norm(transformed_vector_f1_F5) * np.linalg.norm(transformed_vector_f2_F5)))

# Change in length for fibers f1 and f2
change_length_f1_F4 = length_f1_F4 - original_length_f1
change_length_f2_F4 = length_f2_F4 - original_length_f2
change_length_f1_F5 = length_f1_F5 - original_length_f1
change_length_f2_F5 = length_f2_F5 - original_length_f2

# Change in angle between fibers f1 and f2
change_angle_f1_f2_F4 = angle_f1_f2_F4 - angle_f1_f2_original
change_angle_f1_f2_F5 = angle_f1_f2_F5 - angle_f1_f2_original

print('change_length_f1_F4=',change_length_f1_F4)
print('change_length_f2_F4=',change_length_f2_F4)
print('change_length_f1_F5=',change_length_f1_F5)
print('change_length_f2_F5=',change_length_f2_F5)
print('change_angle_f1_f2_F4=',change_angle_f1_f2_F4)
print('change_angle_f1_f2_F5=',change_angle_f1_f2_F5)

