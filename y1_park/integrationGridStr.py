import numpy as np

from integrationVolumeStr import radius

rp = radius[1:]
nradii = len(rp)
z_low = [0.2, 0.35, 0.5]
z_hig = [0.35, 0.5, 0.65]

z_wall = np.repeat(z_low,nradii)
z_wall_numbers = "zo_low  = "+" ".join(map(str, list(z_wall)))
print(z_wall_numbers)

z_wall = np.repeat(z_hig,nradii)
z_wall_numbers = "zo_high = "+" ".join(map(str, list(z_wall)))
print(z_wall_numbers)

rp3 = np.tile(rp,3)
wall_numbers = "radii  = "+" ".join(map(str, list(rp3)))
print(wall_numbers)

wall_numbers = "radius = "+" ".join(map(str,rp))
print(wall_numbers)
