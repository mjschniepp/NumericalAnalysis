'''
 Author: Michael Schniepp
 Date Modified: 11/11/15
 Implementation and Expirimentation of NumPy's FFT package for the purpose of trigonometric polynomial interpolation.
 Test Version # 2
'''

import numpy as np
import matplotlib.pyplot as plt

# The function to generate f_j points
def fx(x):
    y = np.exp( np.sin(x) )
    return y

# The derivative of the function for later comparison
def fxPrime(x):
    y = np.exp( np.sin(x) )*np.cos(x)
    return y
  
# x_j node generation
n = 32
x_nodes = []
for j in range(0,n):
    x_nodes.append(j*2*np.pi / float(n))
    
# f_j points generated
fxj = []
for j in range(0,n):
    fxj.append(fx(x_nodes[j]))  
    
# FFT package implemented and altered to match defenition of c_k
dft = np.fft.fft(fxj)
for i in range(0,n):
    dft[i] = dft[i]/float(n)
    
N = n/2
x = np.arange(0,2*np.pi,0.02)
k = np.arange(1,N)

# Generation of our interpolationg polynomial
p = []
for i in x:
    acc = 0
    for j in k:
        acc += 2*dft[j]*np.exp(1j*j*i) # again note use of 2 times index 0 to N
    acc += float(0.5)*dft[N]*np.cos(N*j*i)
    p.append(acc)
    
# By multiplying ik by our above sum we get the derivative interpolating polynomial
p_prime = []
for i in x:
    acc = 0
    for j in k:
        acc += (1j*j)*(2*dft[j]*np.exp(1j*j*i))
    acc += float(0.5)*dft[N]*np.cos(N*j*i)
    p_prime.append(acc)    
    
# Finalizing P_N            
a0 = dft[0]
p2 = a0 + p
p3 = p_prime

# Error computed subtracting the interpolated polynomial from the real derivative
# when N is doubled the error reduced by a factor of 10,000
z = fxPrime(x)
err = []
for i in range(0,315):
    err.append( z[i] - p3[i] )
 # plotting code
f, axarr = plt.subplots(3, sharex=True)
axarr[0].set_title('P_8')
axarr[0].plot(x_nodes,fxj,'ro',x,p2,x,fx(x))
axarr[1].set_title('P_8 Derivative')
axarr[1].plot(x,p3, x,z)
axarr[2].set_title('Error Plot')
axarr[2].plot(x,err)
plt.show()
