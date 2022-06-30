'''
Fourier Approximations
- Fitting the functions exp(x) and cos(cos(x)) over the interval [0, 2pi) 
  using the Fourier series expansion
- Finding the Fourier coefficients using least-squares method and 
  analyzing their deviation from actual coefficients
- Drawing important inferences through a variety of plots

Mar 11, 2022

By Arjun R
EE20B016
'''

from scipy.integrate import quad
from scipy.linalg import lstsq
import numpy as np
from pylab import *

"""1"""
# Generating the required functions
exp = lambda x: np.exp(x)
coscos = lambda x: np.cos(np.cos(x))
pi = np.pi

# Generating the required plotting ranges
N = 1201
x_ = np.linspace(-2*pi, 4*pi, N)
x_ = x_[:-1]
x = linspace(0, 2*pi, 401)
x = x[:-1] # drop last term to have a proper periodic integral

figure(1)
semilogy(x_, exp(x_), label='True')
expected_exp = np.concatenate((exp(x), exp(x), exp(x)))
semilogy(x_, expected_exp, label='Expected', alpha=0.8)
xlabel(r'$x \rightarrow$')
xticks(x_[::int(N/12)], rotation=90)
ylabel(r'$exp(x) \rightarrow$')
legend(loc = 'upper right')
title('Exponential function')
grid()
show()

figure(2)
plot(x_, coscos(x_), label='True')
expected_coscos = np.concatenate((coscos(x), coscos(x), coscos(x)))
plot(x_, expected_coscos, label='Expected', alpha=0.8)
xlabel(r'$x \rightarrow$')
xticks(x_[::int(N/12)], rotation=90)
ylabel(r'$cos(cos(x)) \rightarrow$')
legend(loc = 'upper right')
title('cos(cos(x)) function')
grid()
show()

"""2"""
# Functions to be integrated
def u(x, f, k):
    return f(x) * cos(k*x)

def v(x, f, k):
    return f(x) * sin(k*x)

"""3""" # Generating the Fourier coefficients using quad()
c_exp = np.zeros(51)
c_exp[0] = quad(u, 0, 2*pi, args=(exp, 0))[0]/(2*pi)
for k in range(1, 26):
    c_exp[2*k-1] = quad(u, 0, 2*pi, args=(exp, k))[0]/pi
    c_exp[2*k] = quad(v, 0, 2*pi, args=(exp, k))[0]/pi

#print(c_exp) 

c_coscos = np.zeros(51)
c_coscos[0] = quad(u, 0, 2*pi, args=(coscos, 0))[0]/(2*pi)
for k in range(1, 26):
    c_coscos[2*k-1] = quad(u, 0, 2*pi, args=(coscos, k))[0]/pi
    c_coscos[2*k] = quad(v, 0, 2*pi, args=(coscos, k))[0]/pi

#print(c_coscos)

figure(3); 
semilogy(np.arange(26), np.abs(c_exp[::2]), 'ro', label=r'$a_n$'); 
semilogy(np.arange(25), np.abs(c_exp[1::2]), 'o', color = 'lightcoral', label=r'$b_n$'); 
xlabel('n'); legend(); title('Semilogy plot for exponential function'); grid(); show()
figure(4); 
loglog(np.arange(26), np.abs(c_exp[::2]), 'ro', label=r'$a_n$'); 
loglog(np.arange(25), np.abs(c_exp[1::2]), 'o', color = 'lightcoral', label=r'$b_n$'); 
xlabel('n'); legend(); title('Loglog plot for exponential function'); grid(); show()
figure(5); 
semilogy(np.arange(26), np.abs(c_coscos[::2]), 'ro', label=r'$a_n$'); 
semilogy(np.arange(25), np.abs(c_coscos[1::2]), 'o', color = 'lightcoral', label=r'$b_n$'); 
xlabel('n'); legend(); title('Semilogy plot for cos(cos(x)) function'); grid(); show()
figure(6); 
loglog(np.arange(26), np.abs(c_coscos[::2]), 'ro', label=r'$a_n$'); 
loglog(np.arange(25), np.abs(c_coscos[1::2]), 'o', color = 'lightcoral', label=r'$b_n$'); 
xlabel('n'); legend(); title('Loglog plot for cos(cos(x)) function'); grid(); show()

"""4"""
#x = linspace(0, 2*pi, 401)
#x = x[:-1]
b1 = exp(x)
A = zeros((400, 51)) # allocate space for A
A[:, 0] = 1 # col 1 is all ones
for k in range(1, 26):
    A[:, 2*k-1] = cos(k*x) # cos(kx) column
    A[:, 2*k] = sin(k*x) # sin(kx) column
    #endfor
c1 = lstsq(A, b1, rcond = None)[0] # the ’[0]’ is to pull out the
# best fit vector. lstsq returns a list.

# Similarly applying least-squares method for cos(cos(x))
b2 = coscos(x)    
c2 = lstsq(A, b2, rcond = None)[0]

#print(c1)
#print(c2)

figure(3) 
semilogy(np.arange(26), np.abs(c1[::2]), 'go', label=r'$est. a_n$'); 
semilogy(np.arange(25), np.abs(c1[1::2]), 'o', color = 'limegreen', label=r'$est. b_n$'); 
"""Including old plots"""
semilogy(np.arange(26), np.abs(c_exp[::2]), 'ro', label=r'$a_n$'); 
semilogy(np.arange(25), np.abs(c_exp[1::2]), 'o', color = 'lightcoral', label=r'$b_n$'); 
xlabel('n'); legend(); title('Semilogy plot for exponential function'); grid(); show()

figure(4)
loglog(np.arange(26), np.abs(c1[::2]), 'go', label=r'$est. a_n$'); 
semilogy(np.arange(25), np.abs(c1[1::2]), 'o', color = 'limegreen', label=r'$est. b_n$'); 
"""Including old plots"""
loglog(np.arange(26), np.abs(c_exp[::2]), 'ro', label=r'$a_n$'); 
semilogy(np.arange(25), np.abs(c_exp[1::2]), 'o', color = 'lightcoral', label=r'$b_n$');  
xlabel('n'); legend(); title('Loglog plot for exponential function'); grid(); show()

figure(5) 
semilogy(np.arange(26), np.abs(c2[::2]), 'go', label=r'$est. a_n$'); 
semilogy(np.arange(25), np.abs(c2[1::2]), 'o', color = 'limegreen', label=r'$est. b_n$'); 
"""Including old plots"""
semilogy(np.arange(26), np.abs(c_coscos[::2]), 'ro', label=r'$a_n$'); 
semilogy(np.arange(25), np.abs(c_coscos[1::2]), 'o', color = 'lightcoral', label=r'$b_n$'); 
xlabel('n'); legend(); title('Semilogy plot for cos(cos(x)) function'); grid(); show()

figure(6)
semilogy(np.arange(26), np.abs(c2[::2]), 'go', label=r'$est. a_n$'); 
semilogy(np.arange(25), np.abs(c2[1::2]), 'o', color = 'limegreen', label=r'$est. b_n$'); 
"""Including old plots"""
semilogy(np.arange(26), np.abs(c_coscos[::2]), 'ro', label=r'$a_n$'); 
semilogy(np.arange(25), np.abs(c_coscos[1::2]), 'o', color = 'lightcoral', label=r'$b_n$'); 
xlabel('n'); legend(); title('Loglog plot for cos(cos(x)) function'); grid(); show()

"""6"""
# Calculating the absolute deviations between both sets of coefficents
c_exp_dev = np.abs(c_exp-c1)
c_coscos_dev = np.abs(c_coscos-c2)
print("Max. absolute deviation between coefficients of exp(x) generated through both methods: ", np.max(c_exp_dev))
print("Max. absolute deviation between coefficients of cos(cos(x)) generated through both methods: ", np.max(c_coscos_dev))

# Plotting the function generated by coefficents estimated through least-squares method
f1 = np.dot(A, c1) # exp(x)
f2 = np.dot(A, c2) # cos(cos(x))
figure(1)
semilogy(x_, exp(x_), label='True')
semilogy(x_, expected_exp, label='Expected')
semilogy(x, f1, 'go', markersize = 2, label = 'Estimate')
xlabel(r'$x \rightarrow$')
xticks(x_[::int(N/12)], rotation=90)
ylabel(r'$exp(x) \rightarrow$')
legend(loc = 'upper right')
title('Exponential function')
grid()
show()

figure(2)
plot(x_, coscos(x_), label='True')
plot(x_, expected_coscos, label='Expected')
plot(x, f2, 'go', markersize = 2, label = 'Estimate')
xlabel(r'$x \rightarrow$')
xticks(x_[::int(N/12)], rotation=90)
ylabel(r'$cos(cos(x)) \rightarrow$')
legend(loc = 'upper right')
title('cos(cos(x)) function')
grid()
show()

