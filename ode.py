import math

def basic_ode(x, t, lam, beta, gamma1, gamma2, gamma3, k_alpha, k_mu, suppression):
    S1 = x[0]
    S2 = x[1]
    I1 = x[2]
    I2 = x[3]
    # ---
    Infections = x[4]

    alpha = -k_alpha * (I2/(I1+I2) - suppression)
    mu = -k_mu * (S2-I1)

    dS1 = beta * (S1+S2+I1+I2) - (mu+gamma1)*S1
    dS2 = mu*S1 - (gamma1+lam)*S2
    dI1 = lam*S2 - (gamma2+alpha)*I1
    dI2 = alpha*I1 - gamma3*I2
    #dR = gamma1*(S1+S2) + gamma2*I1 + gamma3*I2

    dInfections = lam*S2
    dDeaths = gamma1*(S1+S2) + gamma2*I1 + gamma3*I2
    dHIVDeaths = (gamma2-gamma1)*I1 + (gamma3-gamma1)*I2

    return [dS1, dS2, dI1, dI2,       dInfections, dDeaths, dHIVDeaths]
