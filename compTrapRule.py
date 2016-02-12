'''
    This is a script to implement the Composite Trapezoidal Rule Quatrature
    This will approximate the value of a function of a single variable on a
    defined interval. Be mindful of n values used. 
    Author: Michael Schniepp
'''

import numpy as np

def compTrap(a,b,n,func):
    h = (b-a)/float(n)
    nodes = []
    processedNodes = []
    
    for i in range(0,n+1):
        nodes.append(a + i*h)
    for i in range(0,n+1):
        processedNodes.append(func(nodes[i]))
    
    processedNodes[0] = 0.5*processedNodes[0]
    processedNodes[-1] = 0.5*processedNodes[-1]
    
    integralValue = h*sum(processedNodes)
    
    return integralValue

# Test functions:
def f(x):
    y = np.exp(x)*x**2
    return y


    
# Simpsons Rule (Richardsons Extrapolation):
def S(a,b,n,func):
    y = compTrap(a,b,n,func)+(4.0/3.0)*(compTrap(a,b,2*n,func)+compTrap(a,b,n,func))
    return y