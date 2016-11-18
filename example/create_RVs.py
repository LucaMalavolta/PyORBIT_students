import numpy as np

import sys
import os
if os.path.isdir('/Users/malavolta/Astro/CODE/'):
    sys.path.insert(0, '/Users/malavolta/Astro/CODE/pyorbit/')
else:
    sys.path.insert(0, '/home/malavolta/CODE/pyorbit/')

import kepler_exo as kp
import scipy.optimize

P = 33.52000
K = 56.39157243567
f = 3.90901023359
ecoso = 0.5394068701567
esino = -0.2592154221098

e = ecoso**2 + esino**2
o = np.arctan2(esino,ecoso)


G_grav = 6.67398e-11
M_sun = 1.98892e30
M_jup = 1.89813e27
M_ratio = M_sun / M_jup
Mu_sun = 132712440018.9
seconds_in_day = 86400
AU_km = 1.4960 * 10 ** 8
def get_mass(M_star2, M_star1, Period, K1, e0):
    # M_star1, M_star2 in solar masses
    # P in days -> Period is converted in seconds in the routine
    # inclination assumed to be 90 degrees
    # Gravitational constant is given in m^3 kg^-1 s^-2
    # output in m/s
    output = K1 - (2. * np.pi * G_grav * M_sun / 86400.0) ** (1.0 / 3.0) * (1.000 / np.sqrt(1.0 - e0 ** 2.0)) * (
                                                                                                                    Period) ** (
                                                                                                                    -1.0 / 3.0) * (
                      M_star2 * (M_star1 + M_star2) ** (-2.0 / 3.0))
    return output

M_star = 1.00

x0 =1/150.0
M = 333032.217814 * scipy.optimize.fsolve(get_mass, x0, args=(M_star, P, K, e))


print M
print e, o

#[ 265.66443948]
#0.358152406632 -0.447972085322

tref = 6050.0000000000000


x = np.arange(6000,6200,1)
xrand = np.random.normal(x,0.1)
xrand0 = xrand-tref
y = kp.kepler_RV_T0P(xrand0, f, P, K, e, o)
yrand = np.random.normal(y,1.)
tcent = kp.kepler_Tcent_T0P(P, f, e, o)

transit = tref + tcent #+ P*np.arange(-5.0,5.0,1.0)

import matplotlib.pyplot as plt
plt.plot(x,y)
plt.scatter(xrand,yrand)
plt.axvline(transit)
plt.show()


#fileout = open('rv_data_1p.dat','w')
#for xv,yv in zip(xrand,yrand):
#    fileout.write('{0:12f} {1:8f} {2:8f} {3:2d} {3:4d}\n'.format(xv,yv,1.0000, 0, 0))
#fileout.close()

fileout = open('tran_1p.dat','w')
fileout.write('{0:5d} {1:8f} {2:8f}\n'.format(0,transit,0.001))
fileout.close()
