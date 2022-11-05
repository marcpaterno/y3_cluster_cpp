import numpy as np
import hydro_mc

# code to implement Child, Habib, Heitmann et al 2018's eq 18,
# a mass-concentration relation.


# m200 is mass M200m in linear units of solar mass/h^-1
def c_from_m200 (m200, z, omega_m, omega_b, sigma8, h0) :
    # using parameters frrom individual, all row of table 1
    A =  3.44
    b =  430.49
    m = -0.10
    co =  3.19
    # now M_* should be calculated using equations 13-17, but we will approximate
    # it as a linear relation from log(M_*/h^_1) = 12.5 at z=0 and 11.0 at z=1
    Mstar =  12.5 + -1.5*z

    # eq 18 demands M200c, not M200m. So we will use the conversion equations in
    # Ragagnin, Saro, Singh, & Dolag 2021 and code at 
    # https://github.com/aragagnin/hydro_mc/blob/master/examples/sample_mm.py
    # so we will need to import hydro_mc
    # https://github.com/aragagnin/hydro_mc/blob/master/hydro_mc.py
    # Then:
    M_200c = hydro_mc.mass_from_mm_relation(
        '200c','200m', M=m200, a=1./(1+z),
        omega_m= omega_m, omega_b= omega_b, sigma8= sigma8, h0= h0)

    # Now we do eq 18
    mmb = M_200c/Mstar/b
    c200 = A * ( mmb**m * (1+mmb)**(-1*m) -1) + co

    return c200
