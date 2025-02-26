import numpy as np

from integrationVolumeStr import radius

rp = radius#[:1]
nradii = len(rp)
z_low = [0.2, 0.35, 0.5]
z_hig = [0.35, 0.5, 0.65]

l_low = [20, 30, 45, 60]

z_lo, l_lo = np.meshgrid(z_low, l_low, indexing='ij')
z_hi, l_hi = np.meshgrid(z_hig, l_low, indexing='ij')

z_lo_wall = z_lo.flatten() #np.repeat(z_low,nradii)
z_wall_numbers = "zo_low  = "+" ".join(map(str, list(z_lo_wall)))
print(z_wall_numbers)

z_hi_wall = z_hi.flatten() #np.repeat(z_low,nradii)
z_wall_numbers = "zo_high  = "+" ".join(map(str, list(z_hi_wall)))
print(z_wall_numbers)

l_lo_wall = l_lo.flatten() +1 #np.repeat(z_low,nradii)
l_wall_numbers = "lo_bin  = "+" ".join(map(str, list(l_lo_wall)))
print(l_wall_numbers)

# rp3 = np.tile(rp,3)
# wall_numbers = "radii  = "+" ".join(map(str, list(rp3)))
# print(wall_numbers)

# wall_numbers = "radius = "+" ".join(map(str,rp))
# print(wall_numbers)
