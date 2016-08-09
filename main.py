import numpy as np
import matplotlib.pyplot as plt
from ode_rats import rats
from scipy.integrate import odeint

# Initial  conditions
y0 = 2.0

# Time grid for  integration
tt = np.linspace(0, 364, 365)

# Solve the ODE
yy = odeint ( rats, y0, tt )

plt.figure(1)
plt.plot(tt, yy, '.')
plt.xlabel( 'Time [days]' )
plt.ylabel( 'Rats [rats]' )
plt.title( 'ODE_Solution' )
plt.show()
