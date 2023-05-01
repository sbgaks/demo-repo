import numpy as np
import matplotlib.pyplot as plt

# Given parameters
theta4 = lambda theta2: 65 + 0.43*theta2 # Function for output angle
theta2_min = 150 # Minimum input angle
theta2_max = 1650 # Maximum input angle
L1 = 410 # Length of fixed link

# Chebyshev's spacing for three precision points
theta2_cheb = (theta2_min + theta2_max)/2 + (theta2_max - theta2_min)/2*np.cos((np.array([0, 1, 2])+0.5)*np.pi/2)

# Solve for link lengths ratios K1, K2, and K3 using the precision points
A = np.array([[np.cos(np.deg2rad(theta2_cheb[0])), -1, 0],
              [np.cos(np.deg2rad(theta2_cheb[1])), -1, 0],
              [np.cos(np.deg2rad(theta2_cheb[2])), -1, 0]])
b = np.array([L1, 0, 0])
x = np.linalg.solve(A, b)
K1, K2, K3 = x

# Calculate transmission angles for given range of input angles and plot
theta2_range = np.arange(theta2_min, theta2_max+50, 50) # Range of input angles
theta4_range = theta4(theta2_range) # Corresponding output angles
theta3_range = np.arccos((K1**2 + K2**2 - K3**2 - L1**2 - 2*K1*K2*np.cos(np.deg2rad(theta2_range))) / (2*K1*K3)) # Transmission angles
plt.plot(theta2_range, np.rad2deg(theta3_range))
plt.xlabel('Input Angle (deg)')
plt.ylabel('Transmission Angle (deg)')
plt.title('Transmission Angle vs. Input Angle')
plt.show()
