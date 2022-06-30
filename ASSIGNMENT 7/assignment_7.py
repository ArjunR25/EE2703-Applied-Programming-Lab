"""
Circuit Analysis Using Sympy
- Analyses filter circuits using Laplace transforms 
- Analysis involves step response, response to sum of sinusoids, and damped inputs

April 08, 2022

By Arjun R
EE20B016
"""

import sympy as sy
import pylab as p
import scipy.signal as sp

def sympy_to_transferfn(xpr, s=sy.symbols('s')):
    """ Takes in Sympy transfer function polynomial 
        Returns Scipy LTI system expression's Nr and Dr """
    num, den = sy.simplify(xpr).as_numer_denom()
    # simplify() reduces the expression, as_numer_denom() returns a/b as a, b
    p_num_den = sy.Poly(num, s), sy.Poly(den, s)  # numer and denom polynomials
    c_num_den = [p.all_coeffs() for p in p_num_den]  # coefficients of the s-polynomial
    l_num, l_den = [sy.lambdify((), c)() for c in c_num_den]  # convert to floats
    return l_num, l_den

def sympy_to_lti(xpr, s=sy.symbols('s')):
    '''Generates LTI transfer function from sympy expression'''
    N, D = sympy_to_transferfn(xpr)
    return sp.lti(N, D)

def stepresponse(H, t):
    '''Generates step response from system transfer function'''
    N, D = sympy_to_transferfn(H)
    D.append(0) # Multiplying H with 1/s
    H_u = sp.lti(N, D)
    t, y = sp.impulse(H_u, None, t)
    return y

def damped_input(t, w0, decay):
    '''Generates damped sinusoidal input
    
    t - time vector
    w0 - angular frequency
    decay - decay constant'''

    return p.cos(w0*t)*p.exp(-1*decay*t)

def lowpass(R1, R2, C1, C2, G, Vi):
    s = sy.symbols('s')
    A = sy.Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
    b = sy.Matrix([0,0,0,-Vi/R1])
    V = A.inv()*b
    return (A,b,V)

A,b,V = lowpass(10000,10000,1e-9,1e-9,1.586,1)
s = sy.symbols('s')
Vo_LP = V[3]
print(sy.simplify(Vo_LP))
ww = p.logspace(0, 8, 801)
ss = 1j*ww
hf = sy.lambdify(s, Vo_LP, 'numpy')
p.loglog(ww, abs(hf(ss)), lw=2)
p.grid(True)
p.title('Low-Pass Filter Magnitude plot')
p.xlabel(r'$\omega\rightarrow$')
p.show()

def highpass(R1, R3, C1, C2, G, Vi):
    s = sy.symbols('s')
    A = sy.Matrix([[0,1,0,-1/G],[-s*R3*C2/(1+s*R3*C2),0,1,0],[0,G,-G,1],[s*C1+s*C2+1/R1,0,-s*C2,-1/R1]])
    b = sy.Matrix([0,0,0,Vi*s*C1])
    V = A.inv()*b
    return (A,b,V)

A,b,V = highpass(10000,10000,1e-9,1e-9,1.586,1)
s = sy.symbols('s')
Vo_HP = V[3]
print(sy.simplify(Vo_HP))
ww = p.logspace(0, 8, 801)
ss = 1j*ww
hf = sy.lambdify(s, Vo_HP, 'numpy')
p.loglog(ww, abs(hf(ss)), lw=2)
p.grid(True)
p.title('High-Pass Filter Magnitude plot')
p.xlabel(r'$\omega\rightarrow$')
p.show()

# Response to sum of sinusoids

t = p.linspace(0, 1e-3, 10000)
Vi = (p.sin(2000*p.pi*t)+p.cos(2e6*p.pi*t))*(t>0)
p.plot(t, Vi)
p.grid(True)
p.title('Sum of sinusoids')
p.xlabel(r'$t\rightarrow$')
p.show()

LP_Vout_sum = sp.lsim(sympy_to_lti(Vo_LP), Vi, t)[1]
p.plot(t, LP_Vout_sum)
p.grid(True)
p.title('LPF O/P for sum of sinusoids')
p.xlabel(r'$t\rightarrow$')
p.show()

HP_Vout_sum = sp.lsim(sympy_to_lti(Vo_HP), Vi, t)[1]
p.plot(t, HP_Vout_sum)
p.grid(True)
p.title('HPF O/P for sum of sinusoids')
p.xlabel(r'$t\rightarrow$')
p.show()

# Zoomed in
p.plot(t[:100], HP_Vout_sum[:100])
p.grid(True)
p.title('HPF O/P for sum of sinusoids')
p.xlabel(r'$t\rightarrow$')
p.show()

# Step response

t = p.linspace(0, 1e-3, 10000)
Step_LP = stepresponse(Vo_LP, t)
p.plot(t, Step_LP)
p.grid(True)
p.title('Step response of LPF')
p.xlabel(r'$t\rightarrow$')
p.show()

Step_HP = stepresponse(Vo_HP, t)
p.plot(t, Step_HP)
p.grid(True)
p.title('Step response of HPF')
p.xlabel(r'$t\rightarrow$')
p.show()

# Damped HF sinusoid parameters
t = p.linspace(0, 1e-3, 10000)
w0 = 1e7
decay = 5e4
V_damped_HF = damped_input(t, w0, decay)
p.plot(t[:1000], V_damped_HF[:1000])
p.grid(True)
p.title('Damped high frequency sinusoidal input')
p.xlabel(r'$t\rightarrow$')
p.show()

Vout_LP = sp.lsim(sympy_to_lti(Vo_LP), V_damped_HF, t)[1]
p.plot(t[:1000], Vout_LP[:1000])
p.grid(True)
p.title('Damped sinusoid input response of LPF')
p.xlabel(r'$t\rightarrow$')
p.show()

Vout_HP = sp.lsim(sympy_to_lti(Vo_HP), V_damped_HF, t)[1]
p.plot(t[:1000], Vout_HP[:1000])
p.grid(True)
p.title('Damped sinusoid input response of HPF')
p.xlabel(r'$t\rightarrow$')
p.show()

# Damped Low frequency sinusoid
t = p.linspace(0, 1e-1, 10000)
w0 = 1e3
decay = 10
V_damped_LF = damped_input(t, w0, decay)
p.plot(t, V_damped_LF)
p.grid(True)
p.title('Damped low frequency sinusoidal input')
p.xlabel(r'$t\rightarrow$')
p.show()

Vout_LP = sp.lsim(sympy_to_lti(Vo_LP), V_damped_LF, t)[1]
p.plot(t, Vout_LP)
p.grid(True)
p.title('Damped sinusoid input response of LPF')
p.xlabel(r'$t\rightarrow$')
p.show()

Vout_HP = sp.lsim(sympy_to_lti(Vo_HP), V_damped_LF, t)[1]
p.plot(t, Vout_HP)
p.grid(True)
p.title('Damped sinusoid input response of HPF')
p.xlabel(r'$t\rightarrow$')
p.show()

# Zoomed in
p.plot(t[:100], Vout_HP[:100])
p.grid(True)
p.title('Damped sinusoid input response of HPF')
p.xlabel(r'$t\rightarrow$')
p.show()

