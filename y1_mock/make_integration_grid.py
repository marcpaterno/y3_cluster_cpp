import numpy as np

nradii = 2
radius = np.round(np.logspace(np.log10(0.05), np.log10(15), nradii),2)

wall_numbers = " ".join(map(str,radius))
print(wall_numbers)

radius3 = np.tile(radius,3)
wall_numbers = " ".join(map(str, list(radius3)))
print(wall_numbers)

z_wall = np.repeat([0.2, 0.373, 0.51],nradii)
z_wall_numbers = " ".join(map(str, list(z_wall)))
print(z_wall_numbers)

z_wall = np.repeat([0.32, 0.51, 0.64],nradii)
z_wall_numbers = " ".join(map(str, list(z_wall)))
print(z_wall_numbers)
