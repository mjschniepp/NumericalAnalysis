'''
    Solving an ODE via Solving a Tridiagonal Matrix
    Author: Michael Schniepp
    we are solving the following equation:
    u''(x) - pi^2*cos(pi*x)u(x) = pi^2*sin(pi*x)cosh(sin(pi*x))
'''

import numpy as np
import timeit 
import matplotlib.pyplot as plt


# preliminary functions and data

    
def GaussianElim(n, A, b):  
    #Gaussian elimination function   
    for j in range(0,n-1):
        for i in range(j+1,n):
            Mij = A[i,j]/A[j,j]
            for k in range(j+1,n):      #(j+1,n):
                A[i,k] = A[i,k] - Mij*A[j,k]
            b[i] = b[i] - Mij*b[j]
    return(A,b)

def TriGaussianElim(n, A, b):  
    #Gaussian elimination function for tridiagonal  
    # Note: this can be improved and is not optimized
    for j in range(0,n-1):
        for i in range(j+1,n):
            Mij = A[i,j]/A[j,j]
            for k in range(j+1,j+2):
                A[i,k] = A[i,k] - Mij*A[j,k]
            b[i] = b[i] - Mij*b[j]
    return(A,b)


def BackSub(n, A, b):
    # backsubstitution implementation
    x = np.zeros(shape=(len(A),1)) 
    for i in range(n-1,-1,-1):
        x[i]= b[i]
        for j in range(i+1, n):
            x[i] = x[i]-A[i,j]*x[j]
        x[i] = x[i]/A[i,i]
    return x
    
def gaussSolver(n,A,b):
    # solves the system
    GaussianElim(n,A,b)
    x = BackSub(n,A,b)
    return x 

def TriGaussSolver(n,A,b):
    # for tridiagonal matrix
    TriGaussianElim(n,A,b)
    x = BackSub(n,A,b)
    return x   
    
def f(x,n):
    # function f(x) given on the hw
    h = 1./(n + 1)
    y = (-h**2)*(np.pi**2)*(np.sin(np.pi*x))*(np.cosh(np.sin(np.pi*x)))
    return y
    

# generating the matricies to be solved, and the nodes values
def gen_nodes(n):
    h = 1./(n + 1)
    nodes = []
    for i in range(0,n+2):
        nodes.append(i*h)
    return nodes
    
def gen_A(n,nodes):
    # populates our tridiagonal matrix
    A = np.zeros((n,n))
    h = 1./(n + 1)
    for i in range(0,n):
        for j in range(0,n):
            if i == j:
                A[i,j] = -2 - (h**2)*(np.pi**2)*(np.cos(np.pi*nodes[i]))**2
            elif i == j+1:
	        A[i,j] = 1
            elif i == j-1:
                A[i,j] = 1
    return A

def gen_b(n,nodes): 
    # populates our b vector
    b = []    
    for i in range(0,n):
        b.append(f(nodes[i],n))
    return b
    
# various n values and test times:

n = [10,20,40,80,160,320]
time = []
triTime = []
errors = []
h = []

# accuracy and computation time tests
for i in range(0,len(n)):
    startTime = timeit.default_timer()
    nodes = gen_nodes(n[i])
    A = gen_A(n[i],nodes)
    b = gen_b(n[i],nodes)
    g = gaussSolver(n[i],A,b)
    endTime = timeit.default_timer() - startTime
    time.append(endTime)
    
for i in range(0,len(n)):
    startTime = timeit.default_timer()
    nodes = gen_nodes(n[i])
    A = gen_A(n[i],nodes)
    b = gen_b(n[i],nodes)
    g = TriGaussSolver(n[i],A,b)
    endTime = timeit.default_timer() - startTime
    triTime.append(endTime)

for i in range(0,len(n)):
    nodes = gen_nodes(n[i])
    A = gen_A(n[i],nodes)
    b = gen_b(n[i],nodes)
    g = TriGaussSolver(n[i],A,b)
    vals = np.ones(len(nodes))
    for j in range(0,len(g)):
        vals[j] = np.sinh(np.sin(np.pi*nodes[j]))
    err = []
    for k in range(0,len(g)):
        val = abs(vals[k] - g[k])
        err.append(val)
    errors.append(max(err))
    h.append(1./(n[i]+1))
    
    
    

logtn = np.log(time)
logn = np.log(n)
tri_logtn = np.log(triTime)
logEh = np.log(errors)
logh = np.log(h)


f, axarr = plt.subplots(2, sharex=False)
plt.subplots_adjust(hspace=0.5)
axarr[0].set_title('Log-Log of time')
axarr[0].plot(logn,logtn,logn,logtn,'ro')
axarr[0].set_xlabel('log(n)')
axarr[0].set_ylabel('log(tn)')
axarr[0].plot(logn,tri_logtn,logn,tri_logtn,'ro')
axarr[1].set_title('Log-Log of Error')
axarr[1].plot(logh,logEh,logh,logEh,'ro')
axarr[1].set_xlabel('log(h)')
axarr[1].set_ylabel('log(e[h])')
plt.show()
    
scale_gauss = (logtn[-1]-logtn[0]) / (logn[-1] - logn[0])
scale_tri = (tri_logtn[-1]-tri_logtn[0]) / (logn[-1] - logn[0])

print 'Scale of Gaussian:', scale_gauss
print 'Scale of Tridiagonal Solver:', scale_tri

# order of error computations, comes to order 1, but im not sure why
# i think the error of max|u-g| was incorrect
order1 = np.log(errors[0]/errors[1])/np.log(2.)
order2 = np.log(errors[1]/errors[2])/np.log(2.)
order3 = np.log(errors[2]/errors[3])/np.log(2.)
order4 = np.log(errors[3]/errors[4])/np.log(2.)
order5 = np.log(errors[4]/errors[5])/np.log(2.)

from tabulate import tabulate
print tabulate([['10/20',order1],
                ['20/40',order2],
                ['40/80',order3],
                ['80/160',order4],
                ['160/320',order5]], 
                headers=['h','Order'])

    
        