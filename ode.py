import math
import matplotlib.pyplot as plt

def basic_ode(y, t):
    #plt.plot(t, y, ko)
    a = 0.01
    omega = 2.0 * math.pi/365.0
    return a * y * ( 1 + math.sin(omega*t) )
