import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from ode import basic_ode
from scipy.integrate import odeint
from scipy.optimize import fsolve # To solve for initial conditions

firstYear = 2015
nYears = 35

# Initial  conditions
def initial_condition_equations(x, population, prevalence, suppression):
    S1 = x[0]
    S2 = x[1]
    I1 = x[2]
    I2 = x[3]

    e1 = sum(x) - population
    e2 = (I1 + I2) / (S1 + S2 + I1 + I2) - prevalence
    e3 = I2 / (I1 + I2) - suppression
    e4 = S2 - I1

    return (e1**2, e2**2, e3**2, e4**2)

population = 1527000 #2262485
prevalence = 0.222

diagnosed = 0.88
on_treatment = 0.88
vl_suppressed = 0.90
print 'Adults on ART', diagnosed*on_treatment
suppression = diagnosed * on_treatment * vl_suppressed

x0 =  fsolve( initial_condition_equations, (0.85*population, 0,0.15*population,1), args=(population, prevalence, suppression) )
print 'Initial conditions:', x0, ' error is: ', sum(initial_condition_equations(x0, population, prevalence, suppression))

# Time grid for  integration
#T = np.linspace(0, 25, 365) # Units of time are years

# Solve the ODE
lam = 0.095             # infection rate applied to the susceptible and at risk (SDC) population
gamma1 = 1/(65.0-15.0)  # death rate from susceptible population
gamma2 = 1/10.0         # death rate from infected (unsuppressed) population
gamma3 = 1/(60.0-15.0)  # death rate from infected (suppressed) population
beta = 2.15*gamma1      # birth rate per capita

k_alpha = 100   # gain to maintain suppression
k_mu = 0.1      # gain to maintain susceptible and at risk population

x0 = np.append(x0, [0,0,0])        # for new infections counter
xp = x0
X = x0
T = firstYear
new_infections = np.zeros(nYears)
new_deaths = np.zeros(nYears)
new_hiv_deaths = np.zeros(nYears)
for yi,year in enumerate(range(firstYear, firstYear+nYears)):
    Tt = np.linspace(year,year+1,365)
    Xt = odeint(basic_ode, xp, Tt, args=(lam, beta, gamma1, gamma2, gamma3, k_alpha, k_mu, suppression))
    X = np.vstack( (X,Xt) )
    T = np.hstack( (T,Tt) )
    xp = Xt[-1,:]

    # Save and reset new infections counter
    new_infections[yi] = xp[4]
    new_deaths[yi] = xp[5]
    new_hiv_deaths[yi] = xp[6]
    xp[4] = 0
    xp[5] = 0
    xp[6] = 0

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

Pop = np.sum(X[:,0:4], axis=1)
print Pop[-1]
plt.plot(T, Pop)
# Source: https://esa.un.org/unpd/wpp/DataQuery/
plt.plot( range(firstYear, firstYear+nYears+5, 5), [p*1e3 for p in [1527, 1680, 1841, 2024, 2192, 2347, 2482, 2602]], 'ko')

ax.yaxis.set_major_formatter(y_formatter)
plt.xlabel( 'Time [years]' )
plt.ylabel( 'Adult Population' )
plt.tight_layout()

## Figure 3: Constraints
fig = plt.figure(3)
ax = fig.add_subplot(121)

## 3a: fraction suppressed
Pop = np.sum(X, axis=1)
ax.plot(T, X[:,3]/(X[:,2]+X[:,3]))
ax.axhline(y=suppression, color='r')

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

plt.plot(T, (diagnosed+on_treatment) * (X[:,2] + X[:,3]) )

ax.yaxis.set_major_formatter(y_formatter)
plt.xlabel( 'Time [years]' )
plt.ylabel( 'Number on ART' )
plt.tight_layout()


## Figure 5: Prevalence, incidence, new infections, PLHIV
fig = plt.figure(5)

ax = fig.add_subplot(231)
ax.yaxis.set_major_formatter(y_formatter)
ax.plot(T, (X[:,2]+X[:,3]) / (X[:,0]+X[:,1] + X[:,2]+X[:,3]) )
plt.xlabel( 'Time [years]' )
plt.ylabel( 'HIV Prevalence' )

ax = fig.add_subplot(232)
ax.yaxis.set_major_formatter(y_formatter)
ax.plot(T, 100 * lam*X[:,1] / (X[:,0]+X[:,1]) )
plt.xlabel( 'Time [years]' )
plt.ylabel( 'HIV Incidence Rate (%)' )
plt.tight_layout()

ax = fig.add_subplot(233)
ax.yaxis.set_major_formatter(y_formatter)
ax.plot(T, X[:,2]+X[:,3] )
plt.xlabel( 'Time [years]' )
plt.ylabel( 'PLHIV' )

ax = fig.add_subplot(234)
ax.yaxis.set_major_formatter(y_formatter)
ax.plot(new_infections )
plt.xlabel( 'Time [years]' )
plt.ylabel( 'New Infections' )


ax = fig.add_subplot(235)
ax.yaxis.set_major_formatter(y_formatter)
ax.plot(new_deaths )
plt.xlabel( 'Time [years]' )
plt.ylabel( 'New Deaths' )

ax = fig.add_subplot(236)
ax.yaxis.set_major_formatter(y_formatter)
ax.plot(new_hiv_deaths )
plt.xlabel( 'Time [years]' )
plt.ylabel( 'New HIV Deaths' )





plt.show()
