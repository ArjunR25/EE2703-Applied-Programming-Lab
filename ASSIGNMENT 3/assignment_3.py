'''
Fitting Data to Models
-We have to run generate_data.py file to generate a function 
with normally distributed noise added to it and stores it to a file fitting.dat
The program:
-Parses the data file and analyses the data to extract information
-Tries to fit the data to a model using least-squares method
-Studies the effect of noise on the fitting process using different kinds of plots

Feb 18, 2022

By Arjun R
EE20B016
'''
import scipy.special as sp
from scipy.linalg import lstsq
import numpy as np
from pylab import *

data = np.loadtxt('fitting.dat') # Takes in the file data as a 101X10 numpy matrix
t = data[:, 0] # The first column in the file signifies points in time
k = 9 # Number of data columns
N = 101 # Number of rows/data points
A0 = 1.05; B0 = -0.105 # True values for function parameters

"""3&4"""
def g(t, A, B): # Generates the exact function
    y = A*sp.jn(2,t) + B*t
    return y

sigma = np.logspace(-1,-3,9) # returns an numpy array of 9 evenly spaced samples in the logarithmic scale from 1e-1 to 1e-3, both inclusive
plt.figure(figsize = [8, 6])
for i in range(1, k+1):
    plot(data[:, 0], data[:, i], label = r'$\sigma_{}$ = {:.4f}'.format(i, sigma[i-1]))

plot(t, g(t, A0, B0), label = 'true value', color = 'black', linewidth = 3)
xlabel(r'$t \rightarrow$')
ylabel(r'$f(t) + noise \rightarrow$')
legend(loc = 'upper right')
title('Data to be fitted to theory')
show() # Plots the data in all the data columns along with the exact function plot

"""5"""
errorbar(t[::5], data[:, 1][::5], sigma[0], label = r'$Errorbar$', fmt='ro')
plot(t, g(t, A0, B0), label = r'$f(t)$')
xlabel(r'$t \rightarrow$')
legend()
title('Data points for $\sigma = 0.10$ with exact function')
show() # Plots error bar plot for sigma = 0.10 data column

"""6"""
M = np.c_[sp.jn(2, t), t] # Joins two column matrices
p_vec = np.array([A0, B0])
print(np.array_equal(np.dot(M, p_vec).transpose(), g(t, A0, B0))) # Can use np.allclose() if this doesn't return True

"""7"""
A_arr = np.arange(0, 2.1, 0.1)
B_arr = np.arange(-0.2, 0.01, 0.01)
MSE = np.zeros((len(A_arr), len(B_arr)))
for i in range(len(A_arr)):
    for j in range(len(B_arr)):
        MSE[i, j] = (1/N)*sum(np.square(data[:, 1][k] - g(t[k], A_arr[i], B_arr[j])) for k in range(N))

print(A_arr, B_arr, MSE)

"""8"""
X = np.array(A_arr)
Y = np.array(B_arr)
cp = contour(X, Y, MSE, 20,)
clabel(cp, np.linspace(0.025, 0.100, 4))
xlabel(r'$A \rightarrow$')
ylabel(r'$B \rightarrow$')
plot(A0, B0, 'ro')
text(A0, B0, 'Exact Location')
title(r'Contour plot of $\epsilon_{ij}$')
show()

"""9"""
p, resid, rank, sig = lstsq(M, g(t, A0, B0), rcond = None)
print(p)
sigma = np.logspace(-1,-3,9)
error = np.zeros((9, 2))
for col in range(1, 10):
    p, resid, rank, sig = lstsq(M, data[:, col], rcond = None)
    for p_i in range(2):
        error[col-1, p_i] = abs(p[p_i]-p_vec[p_i]) # Calculates absolute value of difference between best estimate and exact value

plot(sigma, error[:, 0], 'ro--', label = 'Aerr', linewidth = 0.5)
plot(sigma, error[:, 1], 'go--', label = 'Berr', linewidth = 0.5)
legend()
xlabel(r'$Noise standard deviation \rightarrow$')
ylabel(r'$Absolute error \rightarrow$')
title('Variation of error with noise')
grid()
show() # Produces a non-linearly varying plot

loglog(sigma, error[:, 0], 'ro', label = 'Aerr')
loglog(sigma, error[:, 1], 'go', label = 'Berr')
errorbar(sigma, error[:, 0], np.std(error[:, 0]), fmt='ro')
errorbar(sigma, error[:, 1], np.std(error[:, 1]), fmt='go')
legend()
xlabel(r'$\sigma_n \rightarrow$')
ylabel(r'$Absolute error \rightarrow$')
title('Variation of error with noise')
grid()
show() # Produces a linearly(approx.) varying plot since sigma was generated with logarithmic spacing