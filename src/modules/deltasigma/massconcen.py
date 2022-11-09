import numpy as np
import hydro_mc
from scipy.interpolate import interp1d

# code to implement Child, Habib, Heitmann et al 2018's eq 18,
# a mass-concentration relation.


# m200 is mass M200m in linear units of solar mass/h^-1
def c_from_m200 (m200, z, omega_m, omega_b, sigma8, h0, mstar,mstar_z) :
    # using parameters from individual, all row of table 1
    A =  3.44
    b =  430.49
    m = -0.10
    co =  3.19
    # row 3, stacked NFW
    #A =  4.61
    #b =  638.65
    #m = -0.07
    #co =  3.59
    # row 4, stacked einasto
    #A =  63.2
    #b =  431.48
    #m = -0.01
    #co =  3.36
    # now M_* should be calculated using equations 13-17, but we will approximate it
    # reproduces Colossus' Child18 model reasonably well
    #Mstarc =  10**(14.76 * .808**z)

    # from Child's eq 13 via mstar.py module
    interp=interp1d(mstar_z, mstar)
    Mstar = interp(z)
    #z = 0.01 # c = 5.02 at 14.0
    #z = 0.33 # c = 4.60 at 14.0
    #z = 0.66 # c = 4.20 at 14.0
    #z = 0.99 # c = 3.90 at 14.0

    # eq 18 demands M200c, not M200m. So we will use the conversion equations in
    # Ragagnin, Saro, Singh, & Dolag 2021 and code at 
    # https://github.com/aragagnin/hydro_mc/blob/master/examples/sample_mm.py
    # so we will need to import hydro_mc
    # https://github.com/aragagnin/hydro_mc/blob/master/hydro_mc.py
    # Then:
    M_200c = hydro_mc.mass_from_mm_relation(
        '200m','200c', M=m200, a=1./(1+z),
        omega_m= omega_m, omega_b= omega_b, sigma8= sigma8, h0= h0)
    # note M_200c should be strictly less then M_200m

    # Now we do eq 18
    mmb = M_200c/Mstar/b
    c200 = A * ( mmb**m * (1+mmb)**(-1*m) -1) + co
    #print(np.log10(m200),np.log10(Mstar),m200/Mstar,c200)

    return c200


# use the mass concentration relation from
    # Ragagnin, Saro, Singh, & Dolag 2021 and code at 
    # so we will need to import hydro_mc
    # https://github.com/aragagnin/hydro_mc/blob/master/hydro_mc.py
    #
# m200 is mass M200m in linear units of solar mass/h^-1
def c_from_m200_ragagnin (m200, z, omega_m, omega_b, sigma8, h0) :

    z= 0.5; m200=0.3e15
    cs = c_from_m200_child (m200, z, omega_m, omega_b, sigma8, h0) 

    c200 = hydro_mc.concentration_from_mc_relation(
            '200m', M=m200, a=1./(1+z), omega_m= omega_m, omega_b= omega_b, sigma8= sigma8, h0= h0)

    print("===================",c200, cs)
    import sys
    sys.tracebacklimit = 0
    raise Exception("here")
    return c200
