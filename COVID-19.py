#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 13:43:26 2018
@author:  (based on a code from Oregon State)
S0 - susceptable
I0 - infected 
R0 - recovered
D - dead 
alpha - the probability of being infected with teh corona viris
https://www.pharmaceutical-technology.com/news/china-wuhan-coronavirus-deaths-12-feb/
beta - the probability of someone recovering the corona virus
https://www.statnews.com/2020/02/04/two-scenarios-if-new-coronavirus-isnt-contained/
delta - the death rate of the corona virus
https://www.worldometers.info/coronavirus/coronavirus-death-rate/

"""
# This loads some pacakges that have arrays etc and the ODE integrators
import scipy, scipy.integrate
import pylab
    

# Parameters
alpha = .0017 
beta = .91
delta = .067

#if you want to include the natural population growths and deaths
#the natual population increase each year
#zeta = .0179
#gamma = .0035


# Initial condition
S0 = 5832710
I0 = 58
R0 = 15

Y0 = [ S0, I0, R0]

tMax = 2

# Time vector for solution
T = scipy.linspace(0, tMax, 10)


# This defines a function that is the right-hand side of the ODEs
# Warning!  Whitespace at the begining of a line is significant!
def rhs(Y, t, alpha, beta, delta):
    '''
    SIR model.
    
    This function gives the right-hand sides of the ODEs.
    '''
    
    # Convert vector to meaningful component vectors
    # Note: Indices start with index 0, not 1!
    S = Y[0]
    I = Y[1]
    R = Y[2]

    # the total number of people in the country of signapore, increasing 
    # by the natural rate of increase
    N = (5832710 * 1.79)
   
    
    # The right-hand sides
    dI =  alpha * I * S  - beta * I - delta * I
    dS = -alpha * I * S 
    dR = beta * I
    
    
    
    # Convert meaningful component vectors into a single vector
    dY = [ dS, dI, dR ]

    return dY




# Integrate the ODE
# Warning!  The ODE solver over-writes the initial value.
# This can be a pain in the ass if you run the ODE multiple times.
# Also, 'args' passes parameters to right-hand-side function.
solution = scipy.integrate.odeint(rhs,
                                  Y0,
                                  T,
                                  args = (alpha, beta, delta))
        
S = solution[:, 0]
I = solution[:, 1]
R = solution[:, 2]

D = delta * I
N = S + I + R 



# Make plots

# Load a plotting package
# PyLab is motivated by Matlab...


# I usually use PyLab for quick plots
# and the Python GnuPlot package for publication

pylab.figure()

pylab.plot(T, S / N,
           T, I / N,
           T, R / N,
           T, D / N)

pylab.xlabel('Time')
pylab.ylabel('Proportion')

pylab.legend([ 'Susceptible', 'Infected', 'Recovered', 'Dead'])

# Actually display the plot
pylab.show()
