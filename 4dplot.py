import numpy as np
import matplotlib.pyplot as plt

# Example 3D data
data = np.random.rand(20, 30, 40)

# Create coordinate grids
z, y, x = np.indices(data.shape)

# Flatten everything so each voxel becomes one point
x = x.flatten()
y = y.flatten()
z = z.flatten()
values = data.flatten()

# Plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

scatter = ax.scatter(
    x, y, z,
    c=values,          # color represents the 4th dimension
    cmap="viridis",    # colormap
    s=10,              # marker size
    alpha=0.7
)

# Colorbar
cbar = fig.colorbar(scatter, ax=ax, pad=0.1)
cbar.set_label("Value")

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.show()
