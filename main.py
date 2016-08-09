import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
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

population = 2262485
prevalence = 0.22
suppression = 0.9 * 0.9 * 0.9
x0 =  fsolve( initial_condition_equations, (0.85*population, 0,0.15*population,1,0), args=(population, prevalence, suppression) )
print 'Initial conditions:', x0, ' error is: ', sum(initial_condition_equations(x0, population, prevalence, suppression))

# Time grid for  integration
#T = np.linspace(0, 25, 365) # Units of time are years

# Solve the ODE
lam = 0.073         # infection rate applied to the susceptible and at risk (SDC) population
gamma1 = 1/65.0     # death rate from susceptible population
gamma2 = 1/10.0     # death rate from infected (unsuppressed) population
gamma3 = 1/60.0     # death rate from infected (suppressed) population
beta = 2.5*gamma1     # birth rate per capita

k_alpha = 100       # gain to maintain suppression
k_mu = 0.1          # gain to maintain susceptible and at risk population

x0 = np.append(x0, 0)        # for new infections counter
xp = x0
X = x0
T = 0
new_infections = np.zeros(25)
for year in range(25):
    Tt = np.linspace(year,year+1,365)
    Xt = odeint(basic_ode, xp, Tt, args=(lam, beta, gamma1, gamma2, gamma3, k_alpha, k_mu, suppression))
    X = np.vstack( (X,Xt) )
    T = np.hstack( (T,Tt) )
    xp = Xt[-1,:]

    # Save and reset new infections counter
    new_infections[year] = xp[5]
    xp[5] = 0

#X = odeint(basic_ode, x0, T, args=(lam, beta, gamma1, gamma2, gamma3, k_alpha, k_mu, suppression))

## Figure 1
y_formatter = ScalarFormatter(useOffset=False)

fig = plt.figure(1)
ax = fig.add_subplot(111)

labels = ['Susceptible (Not at risk)', 'Susceptible (at risk)', 'Infected (unsuppressed)', 'Infected (suppressed)', 'Died']

for Xi, label in zip(X.T, labels):
    ax.plot(T, Xi, label=label)

ax.yaxis.set_major_formatter(y_formatter)
plt.xlabel( 'Time [years]' )
plt.ylabel( 'Population' )
plt.title( 'ODESolution' )
plt.legend()
plt.tight_layout()

## Figure 2: Population
fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(y_formatter)

Pop = np.sum(X[:,:-1], axis=1)
print Pop[-1]
plt.plot(T, Pop)

ax.yaxis.set_major_formatter(y_formatter)
plt.xlabel( 'Time [years]' )
plt.ylabel( 'Population' )
plt.tight_layout()

## Figure 3: Constraints
fig = plt.figure(3)
ax = fig.add_subplot(121)

## 3a: fraction suppressed
Pop = np.sum(X, axis=1)
ax.plot(T, X[:,3]/(X[:,2]+X[:,3]))
ax.axhline(y=suppression, xmin=min(T), xmax=max(T), color='r')

ax.yaxis.set_major_formatter(y_formatter)
plt.xlabel( 'Time [years]' )
plt.ylabel( 'Fraction suppressed amongst infected' )

## 3a: size of susceptible & at-risk population
ax = fig.add_subplot(122)

Pop = np.sum(X, axis=1)
ax.plot(T, X[:,1], label='Susceptible (at risk)')
ax.plot(T, X[:,2], label='Infected (unsuppressed)')

ax.yaxis.set_major_formatter(y_formatter)
plt.xlabel( 'Time [years]' )
plt.ylabel( 'At risk population constraint' )
plt.tight_layout()


## Figure 4: Suppressed population
fig = plt.figure(4)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(y_formatter)

plt.plot(T, X[:,3])

ax.yaxis.set_major_formatter(y_formatter)
plt.xlabel( 'Time [years]' )
plt.ylabel( 'Suppressed Population' )
plt.tight_layout()


## Figure 5: Prevalence, incidence, new infections, PLHIV
fig = plt.figure(5)

ax = fig.add_subplot(141)
ax.yaxis.set_major_formatter(y_formatter)
ax.plot(T, (X[:,2]+X[:,3]) / (X[:,0]+X[:,1] + X[:,2]+X[:,3]) )
plt.xlabel( 'Time [years]' )
plt.ylabel( 'HIV Prevalence' )

ax = fig.add_subplot(142)
ax.yaxis.set_major_formatter(y_formatter)
ax.plot(T, lam*X[:,1] / (X[:,0]+X[:,1]) )
plt.xlabel( 'Time [years]' )
plt.ylabel( 'HIV Incidence' )
plt.tight_layout()

ax = fig.add_subplot(143)
ax.yaxis.set_major_formatter(y_formatter)
ax.plot(new_infections )
plt.xlabel( 'Time [yers]' )
plt.ylabel( 'New Infections' )

ax = fig.add_subplot(144)
ax.yaxis.set_major_formatter(y_formatter)
ax.plot(T, X[:,2]+X[:,3] )
plt.xlabel( 'Time [yers]' )
plt.ylabel( 'PLHIV' )





plt.show()
