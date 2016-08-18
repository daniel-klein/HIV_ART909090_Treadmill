# Data points from http://www.avert.org/professionals/hiv-around-world/sub-saharan-africa/botswana

import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from ode_cd4 import cd4_progression
from scipy.integrate import odeint
from scipy.optimize import fsolve # To solve for initial conditions

firstYear = 2015
nYears = 500

x0 =  [0.6, 0.4, 0, 0, 0, 0, 0, 0]

# Solve the ODE
T = np.linspace(firstYear, firstYear+nYears, 12*nYears)
X = odeint(cd4_progression, x0, T) # , args=(lam, beta, gamma1, gamma2, gamma3, k_alpha, k_mu, suppression)

print 'ODE Done'

cum = X[:,7]
s = 0
N = 100000
for trial in range(N):
    r = random.random()
    idx = next(i for (i,c) in enumerate(cum) if r<c)
    s += T[idx]-T[0]

print 'Mean survival is %f years' % (s/N)

# Try with sum and trapz, just for fun ... should match MC :/
print np.sum( [ (t-T[0])*p for (t,p) in zip(T[:-1], np.diff(cum)) ] )
print np.trapz([ (t-T[0])*p for (t,p) in zip(T[:-1], np.diff(cum)) ])

## Figure 1

y_formatter = ScalarFormatter(useOffset=False)

fig = plt.figure(1)
ax = fig.add_subplot(111)

labels = ['>500', '350-500', '250-349', '200-249', '100-199', '50-99', '<50', 'Died']


for Xi, label in zip(X.T, labels):
    ax.plot(T, Xi, label=label)

ax.yaxis.set_major_formatter(y_formatter)
plt.xlim([firstYear, firstYear+25])
plt.ylim([0, 1])
plt.xlabel( 'Time [years]' )
plt.ylabel( 'Population' )
plt.title( 'ODESolution' )
plt.legend()
plt.tight_layout()

plt.show()
