import numpy as np

## set the grid for Rp values
r0, rend, nradii = 0.05, 5, 10
rp = np.logspace(np.log10(r0), np.log10(rend), nradii)
# rp = np.append(0.001, rp)
radius = np.round(np.append(np.array([0.001]), rp),3)

l_low = [20.0, 30.0, 45.0, 60.0]
l_hig = [30.0, 45.0, 60.0, 250]

# lt_low = [5.0]*len(l_low)
# lt_hig = [500]*len(l_low)

lt_low = [int(l/4.) for l in l_low]
lt_hig = [int(l*4.) for l in l_hig]

z_low = [0.15]*len(l_low)
z_hig = [0.75]*len(l_low)

m_low = [28.]*len(l_low)
m_hig = [37.]*len(l_low)

r_low = radius[:-1]
r_hig = radius[1:]

nlam = len(l_low)

def make_integration_volume_str(xlow, xhig, name, nr=nradii):
    xlu, xlb = [], []
    for xl, xh in zip(xlow, xhig):
        xlu += nr * [xl]
        xlb += nr * [xh]
    
    xwall = name+"_low  = "+"    ".join(map(str, xlu))
    xwall+="\n"+name+"_high = "+"    ".join(map(str,xlb))+"\n"
    print(xwall)

def make_radii_tiles(xlow, xhig, name="rp"):
    xlu = np.tile(xlow, nlam)
    xlb = np.tile(xhig, nlam)
    xwall = name+"_low  = "+"    ".join(map(str, xlu))
    xwall+="\n"+name+"_high = "+"    ".join(map(str,xlb))+"\n"
    print(xwall)

def generate_pattern():
    make_integration_volume_str(l_low, l_hig, "lo")
    make_integration_volume_str(lt_low, lt_hig, "lt")
    make_integration_volume_str(z_low, z_hig, "zt")
    make_integration_volume_str(m_low, m_hig, "lnm")
    make_radii_tiles(r_low, r_hig)

if __name__ == "__main__":
    generate_pattern()
