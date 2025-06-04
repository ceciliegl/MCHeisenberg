import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
import matplotlib as mpl
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import equconfig

def color_function(x, y, z):
    # Calculate hue from xy components and brightness from z component
    xy_angle = np.arctan2(y, x)
    hue = (xy_angle + np.pi) / (2 * np.pi)  # Normalize to [0, 1]
    
    if z >= 0:
        # Interpolate linearly between the HSV color and white
        rgb = hsv_to_rgb([hue, 1, 1])
        color = rgb * (1 - z) + np.array([1, 1, 1]) * z
    else:
        # Interpolate linearly between the HSV color and black
        rgb = hsv_to_rgb([hue, 1, 1])
        color = rgb * (1 + z)
        
    return color

# Parameterize the sphere
phi = np.linspace(0, np.pi, 100)
theta = np.linspace(0, 2 * np.pi, 100)
phi, theta = np.meshgrid(phi, theta)

x = np.sin(phi) * np.cos(theta)
y = np.sin(phi) * np.sin(theta)
z = np.cos(phi)

# Generate colors for each point
colors = np.zeros((x.shape[0], x.shape[1], 3))
for i in range(x.shape[0]):
    for j in range(x.shape[1]):
        colors[i, j] = color_function(x[i, j], y[i, j], z[i, j])

# Plot the sphere
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# We have to flatten the arrays for plotting with scatter
ax.plot_surface(x, y, z, facecolors=colors, rstride=1, cstride=1, antialiased=True)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')




# Create a new figure
plt.figure(figsize=(8, 8))

x = np.linspace(-1, 1, 100)
y = np.linspace(-1, 1, 100)
X, Y = np.meshgrid(x, y)

theta = np.arctan2(Y, X).flatten()

print(len(theta))

colors = equconfig.colors(np.cos(2*theta), np.sin(2*theta), np.zeros(len(theta)))
print(colors.shape)

plt.imshow(colors , aspect='equal', extent=(-1.5,1.5,-1.5,1.5))


# Show the plot
plt.show()



