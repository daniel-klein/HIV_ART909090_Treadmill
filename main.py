import numpy as np
import matplotlib.pyplot as plt
from ode import basic_ode
from scipy.integrate import odeint
from scipy.optimize import fsolve # To solve for initial conditions

# Initial  conditions
def initial_condition_equations(x, population, prevalence, suppression):
    S1 = x[0]
    S2 = x[1]
    I1 = x[2]
    I2 = x[3]
    R  = x[4]

    e1 = sum(x) - population
    e2 = R
    e3 = (I1 + I2) / (S1 + S2 + I1 + I2) - prevalence
    e4 = I2 / (I1 + I2) - suppression
    e5 = S2 - I1

    return (e1**2, e2**2, e3**2, e4**2, e5**2)

population = 1e5
prevalence = 0.15
suppression = .9**3
x0 =  fsolve( initial_condition_equations, (0.85*population, 0,0.15*population,1,0), args=(population, prevalence, suppression) )
print 'Initial conditions:', x0, ' error is: ', sum(initial_condition_equations(x0, population, prevalence, suppression))

# Time grid for  integration
T = np.linspace(0, 1, 365) # Units of time are years

# Solve the ODE
lam = 0.1            # infection rate applied to the susceptible and at risk (SDC) population
beta = 1e-3          # birth rate per capita
gamma1 = 1/65.0     # death rate from susceptible population
gamma2 = 1/10.0     # death rate from infected (unsuppressed) population
gamma3 = 1/60.0     # death rate from infected (suppressed) population

k_alpha = 0.1
k_mu = 0.1

X = odeint(basic_ode, x0, T, args=(lam, beta, gamma1, gamma2, gamma3, k_alpha, k_mu, suppression))

print typeof(X)

plt.figure(1)
plt.plot(T, X, '.')
plt.xlabel( 'Time [years]' )
plt.ylabel( 'Population' )
plt.title( 'ODESolution' )
plt.show()
