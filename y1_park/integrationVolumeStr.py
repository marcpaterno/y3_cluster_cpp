import numpy as np

## set the grid for Rp values
r0, rend, nradii = 0.05, 5., 10
rp = np.logspace(np.log10(r0), np.log10(rend), nradii)
# rp = np.append(0.001, rp)
radius = np.round(np.log(np.append(1e-5, rp)), 3)
# 

l_low = [20.0, 30.0, 45.0, 60.0]
l_hig = [30.0, 45.0, 60.0, 200]

# lt_low = [5.0]*len(l_low)
# lt_hig = [500]*len(l_low)

lt_low = [7.76, 12.0, 18.59, 20.76]
lt_hig = [38.0, 53.0, 71.0, 178.0]

# lt_hig = np.where(lt_hig > 200., 200, lt_hig)

# lt_low = [5.0]
# lt_hig = [200]

z_low = [0.15]*len(lt_low)
z_hig = [0.75]*len(lt_low)

m_low = [28.]*len(lt_low)
m_hig = [37.]*len(lt_low)

r_low = radius[:-1]
r_hig = radius[1:]

nlam = len(l_low)
# nradii = 1
def make_integration_volume_str(xlow, xhig, name, nr=nradii):
    xlu, xlb = [], []
    for xl, xh in zip(xlow, xhig):
        xlu += nr * [np.round(xl,3)]
        xlb += nr * [np.round(xh,3)]
    
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
    make_integration_volume_str(np.log(l_low), np.log(l_hig), "logLo")
    make_integration_volume_str(np.log(lt_low), np.log(lt_hig), "logLt")
    make_integration_volume_str(z_low, z_hig, "zt")
    make_integration_volume_str(m_low, m_hig, "lnm")
    make_radii_tiles(r_low, r_hig)

if __name__ == "__main__":
    generate_pattern()
