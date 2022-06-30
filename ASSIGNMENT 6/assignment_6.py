"""
The Laplace Transform
- Analyses LTI systems using scipy's signal processing library
- Analyses Forced Damped Oscillators, Coupled Oscillators and Filter Circuits

Mar 27, 2022

By Arjun R
EE20B016
"""

from pylab import *
import scipy.signal as sp


# Generates transfer function X(s) after solving oscillator equation subject to initial conditions
def transfer_fn(freq, decay): 
    """
    freq - Driving frequency
    decay - Decay constant
    """
    den = polymul([1, 0, 2.25], [1, 2*decay, decay**2 + freq**2])
    return sp.lti([1, decay], den)

'''1'''
t, x = sp.impulse(transfer_fn(freq=1.5, decay=0.5), None, linspace(0, 50, 501)) 
# Calculates impulse response of the transfer function
plot(t, x)
xlabel(r'$t\rightarrow$')
ylabel(r'$x(t)\rightarrow$')
grid()
title('Forced Damped Oscillator with decay = 0.5')
show()

'''2'''
t, x = sp.impulse(transfer_fn(freq=1.5, decay=0.05), None, linspace(0, 50, 501))
plot(t, x)
xlabel(r'$t\rightarrow$')
ylabel(r'$x(t)\rightarrow$')
grid()
title('Forced Damped Oscillator with decay = 0.05')
show()


'''3'''
def f(t, freq, decay): # Forcing function
    return cos(freq*t)*exp(-1*decay*t)*(t>0)

def sys_transfer_fn():
    return sp.lti([1], [1, 0, 2.25])

t = linspace(0, 200, 2001)
fig, ax = subplots(3, 2, figsize = (15, 20))
fig.suptitle('Forced Damped Oscillators (decay = 0.05)', y = 0.9, fontsize = 15)
freq_range = arange(1.4, 1.65, 0.05) # Increasing driving frequencies in steps of 0.05
color = iter(cm.plasma(linspace(0, 1, 5)))
for i in range(5):
    t, x, svec = sp.lsim(sys_transfer_fn(), f(t, freq_range[i], 0.05), t)
    ax[i//2, i%2].plot(t, x, c=next(color))
    ax[i//2, i%2].set_xlabel(r'$t(s)$')
    ax[i//2, i%2].set_ylabel(r'$x$')
    ax[i//2, i%2].set_title(f'freq. = {freq_range[i]} rad/s')
    ax[i//2, i%2].grid()

fig.delaxes(ax[2, 1])
show()

'''4'''
X = sp.lti([1, 0, 2], [1, 0, 3, 0])
Y = sp.lti([2], [1, 0, 3, 0])
t = linspace(0, 20, 201)
t, x = sp.impulse(X, None, t)
t, y = sp.impulse(Y, None, t)

plot(t, x, label='x(t)')
plot(t, y, label='y(t)')
xlabel(r'$t\rightarrow$')
title("Coupled Oscillations")
legend()
grid()
show()

'''5'''
H_ckt = sp.lti([1], [1e-12, 1e-4, 1]) # Transfer function for LP filter circuit
w, mag, phi = H_ckt.bode()
figure(figsize=(12, 5))
subplot(1, 2, 1); semilogx(w, mag)
xlabel(r'$\omega\rightarrow$')
ylabel(r'$20log|H(j\omega)|\rightarrow$')
title('Bode Magnitude Plot'); grid()

subplot(1, 2, 2); semilogx(w, phi)
xlabel(r'$\omega\rightarrow$')
ylabel(r'$\angle H(j\omega)$ (in deg)$\rightarrow$')
title('Bode Phase Plot'); grid(); show()

'''6'''
t = linspace(0.5e-6, 3e-2, 60001)
vi = (cos(1e3*t) - cos(1e6*t)) * (t>0)
t, vo, _ = sp.lsim(H_ckt, vi, t)

plot(t[:61], vo[:61], 'r')
xlabel('t')
ylabel(r'$v_o(t)$')
title('Short term response')
grid()
show()
plot(t, vo, 'r')
xlabel('t')
ylabel(r'$v_o(t)$')
title('Long term response')
grid()
show()

T = linspace(0, 0.030, 1000)
plot(T, cos(1000*T), 'r')
xlabel('t')
title(r'cos(10^3t)')
grid()
show()