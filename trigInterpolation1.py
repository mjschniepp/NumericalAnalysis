'''
 Author: Michael Schniepp
 Date Modified: 11/11/15
 Implementation and Expirimentation of NumPy's FFT package
 Test Version # 2
 Problem Number 5
'''

import numpy as np
import matplotlib.pyplot as plt

# initial data array
A = np.array( [6.000000000000000,10.242640687119284,2.000000000000000,
              -2.585786437626905,2.000000000000000,1.757359312880716,
              -6.000000000000000,-5.414213562373098])
      
# generating x_j nodes   
n = len(A)                   
x_nodes = []
for i in range(0,n):
    x_nodes.append(i*(2.0*np.pi/float(n)))              
        
# implementation of FFT , altered by dividing by n to match defenition 
# of c_k      
y = np.fft.fft(A)
for i in range(0,len(y)):
    y[i] = y[i]/float(n)
     
 
N = n/2
a0 = y[0] # a_0 term 
x = np.arange(0,2*np.pi,0.02)
k = np.arange(1,N)

# Here we compute the sum and last cos component of P_N
p = []
for i in x:
    acc = 0
    for j in k:
        acc += 2*y[j]*np.exp(1j*j*i)     # Note instead of summing -N to N
    acc += float(0.5)*y[N]*np.cos(N*j*i) # we do 2* the sum from 0 to N
    p.append(acc)
         
# Finalizing P_N by adding a_0
p2 = a0 + p   
                   
# Plot to Visualize
plt.title('P_8 Interpolating f_j points')                                                                                                                                                                                                                                    
plt.plot(x,p2,x_nodes,A,'ro')
plt.show()              
              
