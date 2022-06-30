"""
The Digital Fourier Transform
- Computes DFTs for various functions
- Plots their magnitude and phase spectra plots

May 12, 2022

By Arjun R
EE20B016
"""

from pylab import *

# Checking accuracy of DFT
x = rand(100)
X = fft(x)
y = ifft(X)
#print(c_[x,y])
print("Error between fft and ifft:",abs(x-y).max())

# Spectrum of sin(5t)
N = 128
x = linspace(0, 2*pi, N)
y = sin(5*x)
Y = fft(y)
figure()
subplot(211); plot(abs(Y), lw=2)
ylabel(r'$|Y|$', size=16)
title(r'Spectrum of $sin(5t)$')
grid(True)
subplot(212); plot(unwrap(angle(Y)), lw=2)
ylabel(r'Phase of $Y$', size=16)
xlabel(r'$k$', size=16)
grid(True)
show()

N = 128
x = linspace(0, 2*pi, N+1); x = x[:-1] # Removing 2pi to avoid ambiguity
y = sin(5*x)
Y = fftshift(fft(y))/N # Using fftshift and normalizing
w = linspace(-N/2, N/2 -1, N)
figure()
subplot(211); plot(w, abs(Y), lw=2)
xlim([-10, 10])
ylabel(r'$|Y|$', size=16)
title(r'Spectrum of $sin(5t)$')
grid(True)
subplot(212); plot(w, angle(Y), 'ro', lw=2)
ii = where(abs(Y)>1e-3) 
# To highlight the points where spectrum magnitude is significant
plot(w[ii], angle(Y[ii]), 'go', lw=2)
xlim([-10, 10])
ylabel(r'Phase of $Y$', size=16)
xlabel(r'$k$', size=16)
grid(True)
show()

# Spectrum of (1 + 0.1cost)cos10t
t=linspace(0,2*pi,129);t=t[:-1]
y=(1+0.1*cos(t))*cos(10*t)
Y=fftshift(fft(y))/128.0
w=linspace(-64,63,128)
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-15,15])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii = where(abs(Y)>1e-3)
plot(w[ii], angle(Y[ii]), 'go', lw=2)
xlim([-15,15])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()

# increasing number of samples to get fine peaks
t=linspace(-4*pi,4*pi,513);t=t[:-1]
y=(1+0.1*cos(t))*cos(10*t)
Y=fftshift(fft(y))/512.0
w=linspace(-64,64, 513); w = w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-15,15])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii = where(abs(Y)>1e-3)
plot(w[ii], angle(Y[ii]), 'go', lw=2)
xlim([-15,15])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()

# Qn. 2 - Cubic sinusoids
t=linspace(0,2*pi,513);t=t[:-1]
y1 = sin(t)**3
y2 = cos(t)**3
Y1 = fftshift(fft(y1))/512.0
Y2 = fftshift(fft(y2))/512.0
w = linspace(-64, 64, 513); w = w[:-1]
figure()
subplot(221); plot(w, abs(Y1), lw=2)
title(r'Spectrum of $sin^3(t)$')
ylabel(r'$|Y|$', size=16)
xlim([-5, 5])
grid(True)
subplot(222); plot(w, abs(Y2), lw=2)
title(r'Spectrum of $cos^3(t)$')
xlim([-5, 5])
grid(True)
subplot(223); plot(w, angle(Y1), 'ro', lw=2)
ii1 = where(abs(Y1)>1e-3)
plot(w[ii1], angle(Y1[ii1]), 'go', lw=2)
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
xlim([-5, 5])
grid(True)
subplot(224); plot(w, angle(Y2), 'ro', lw=2)
ii2 = where(abs(Y2)>1e-3)
plot(w[ii2], angle(Y2[ii2]), 'go', lw=2)
xlabel(r"$\omega$",size=16)
xlim([-5, 5])
grid(True)
show()

# Qn. 3 - FM signal
t=linspace(-4*pi,4*pi,513);t=t[:-1]
y=cos(20*t+5*cos(t)) 
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513); w = w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-35,35])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of FM signal")
grid(True)
ii = where(abs(Y)>1e-3)
subplot(2,1,2)
plot(w[ii],angle(Y[ii]),'ro',lw=2)
xlim([-35,35])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()

# Qn. 4 - Gaussian function

def Gaussian_exp(x):
    """Computes Gaussian function"""
    return exp(-x**2/2)

def Gauss_fft(w):
    """Computes fft of Gaussian"""
    return (1/sqrt(2*pi))*exp(-w**2/2)

def dft_acc(tolerance = 1e-6, N = 128, T = 8*pi):
    """Computes DFT of Gaussian upto specified value of tolerance"""
    err = tolerance + 1 # to enter the loop
    Y_old = 0

    while err > tolerance:
        x = linspace(-T/2, T/2, N+1)[:-1]
        w = linspace(-pi*N/T, pi*N/T, N+1)[:-1]
        y = Gaussian_exp(x)
        Y = fftshift(fft(ifftshift(y)))*T/(2*pi*N)
        err = max(abs(Y[::2]-Y_old))
        Y_old = Y
        T *= 2
        N *= 2
    
    true_final_error = max(abs(Y - Gauss_fft(w)))

    print("True error: ", true_final_error)
    print("samples = ",N//2, ", time period = pi*", (T/pi)//2)

    figure()
    subplot(2,1,1)
    plot(w,abs(Y),lw=2)
    xlim([-5,5])
    ylabel(r"$|Y|$",size=16)
    title(r"Estimated fft of Gaussian")
    grid(True)
    ii = where(abs(Y)>1e-3)
    subplot(2,1,2)
    plot(w,angle(Y),'ro',lw=2)
    plot(w[ii],angle(Y[ii]),'go',lw=2)
    xlim([-5,5])
    ylim([-2,2])
    ylabel(r"Phase of $Y$",size=16)
    xlabel(r"$\omega$",size=16)
    grid(True)
    show()

    # True fft
    Y1 = Gauss_fft(w)
    figure()
    subplot(2,1,1)
    plot(w,abs(Y1),lw=2)
    xlim([-5,5])
    ylabel(r"$|Y|$",size=16)
    title(r"True fft of Gaussian")
    grid(True)
    ii = where(abs(Y1)>1e-3)
    subplot(2,1,2)
    plot(w,angle(Y1),'ro',lw=2)
    plot(w[ii],angle(Y1[ii]),'go',lw=2)
    xlim([-5,5])
    ylim([-2,2])
    ylabel(r"Phase of $Y$",size=16)
    xlabel(r"$\omega$",size=16)
    grid(True)
    show()

dft_acc()

