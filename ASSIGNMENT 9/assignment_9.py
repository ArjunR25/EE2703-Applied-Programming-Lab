"""
Spectra of Non-periodic Signals
- Computes spectra for non-periodic signals using Hamming window
- Observes time variation of localized DFT spectra using surface plots

May 12, 2022

By Arjun R
EE20B016
"""

from pylab import *
import random
import mpl_toolkits.mplot3d.axes3d as p3

"""Example 1"""
t=linspace(-pi,pi,65);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
y=sin(sqrt(2)*t)
y[0]=0 # the sample corresponding to -tmax should be set zero
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/64.0
w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
#savefig("fig10-1.png")
show()

"""Example 2"""
t1=linspace(-pi,pi,65);t1=t1[:-1]
t2=linspace(-3*pi,-pi,65);t2=t2[:-1]
t3=linspace(pi,3*pi,65);t3=t3[:-1]
# y=sin(sqrt(2)*t)
figure(2)
plot(t1,sin(sqrt(2)*t1),'b',lw=2)
plot(t2,sin(sqrt(2)*t2),'r',lw=2)
plot(t3,sin(sqrt(2)*t3),'r',lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)$")
grid(True)
#savefig("fig10-2.png")
show()

"""Example 3"""
# sin(root(2)t) with t wrapping every 2pi
t1=linspace(-pi,pi,65);t1=t1[:-1]
t2=linspace(-3*pi,-pi,65);t2=t2[:-1]
t3=linspace(pi,3*pi,65);t3=t3[:-1]
y=sin(sqrt(2)*t1)
figure(3)
plot(t1,y,'bo',lw=2)
plot(t2,y,'ro',lw=2)
plot(t3,y,'ro',lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)$ with $t$ wrapping every $2\pi$ ")
grid(True)
#savefig("fig10-3.png")
show()

"""Example 4"""
# Digital ramp spectrum
t=linspace(-pi,pi,65);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
y=t
y[0]=0 # the sample corresponding to -tmax should be set zero
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/64.0
w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
figure()
semilogx(abs(w),20*log10(abs(Y)),lw=2)
xlim([1,10])
ylim([-20,0])
xticks([1,2,5,10],["1","2","5","10"],size=16)
ylabel(r"$|Y|$ (dB)",size=16)
title(r"Spectrum of a digital ramp")
xlabel(r"$\omega$",size=16)
grid(True)
#savefig("fig10-4.png")
show()

"""Example 5"""
# sin(root(2)t) with t wrapping every 2pi after applying Hamming window
t1=linspace(-pi,pi,65);t1=t1[:-1]
t2=linspace(-3*pi,-pi,65);t2=t2[:-1]
t3=linspace(pi,3*pi,65);t3=t3[:-1]
n=arange(64)
wnd=fftshift(0.54+0.46*cos(2*pi*n/63))
y=sin(sqrt(2)*t1)*wnd
figure(3)
plot(t1,y,'bo',lw=2)
plot(t2,y,'ro',lw=2)
plot(t3,y,'ro',lw=2)
ylabel(r"$y$",size=16)
xlabel(r"$t$",size=16)
title(r"$\sin\left(\sqrt{2}t\right)\times w(t)$ with $t$ wrapping every $2\pi$ ")
grid(True)
#savefig("fig10-5.png")
show()

"""Example 6"""
# Using Hamming window for DFT computation
t=linspace(-pi,pi,65);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
n=arange(64)
wnd=fftshift(0.54+0.46*cos(2*pi*n/63))
y=sin(sqrt(2)*t)*wnd
y[0]=0 # the sample corresponding to -tmax should be set zero
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/64.0
w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-8,8])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-8,8])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
#savefig("fig10-6.png")
show()

"""Example 7"""
# Increasing sampling
t=linspace(-4*pi,4*pi,257);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
n=arange(256)
wnd=fftshift(0.54+0.46*cos(2*pi*n/256))
y=sin(sqrt(2)*t)
# y=sin(1.25*t)
y=y*wnd
y[0]=0 # the sample corresponding to -tmax should be set zero
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/256.0
w=linspace(-pi*fmax,pi*fmax,257);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),'b',w,abs(Y),'bo',lw=2)
xlim([-4,4])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y$",size=16)
12
xlabel(r"$\omega$",size=16)
grid(True)
#savefig("fig10-7.png")
show()

# Trying for sin(1.25t)
t=linspace(-4*pi,4*pi,257);t=t[:-1]
dt=t[1]-t[0];fmax=1/dt
n=arange(256)
wnd=fftshift(0.54+0.46*cos(2*pi*n/256))
y=sin(1.25*t)
# y=sin(1.25*t)
y=y*wnd
y[0]=0 # the sample corresponding to -tmax should be set zeroo
y=fftshift(y) # make y start with y(t=0)
Y=fftshift(fft(y))/256.0
w=linspace(-pi*fmax,pi*fmax,257);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),'b',w,abs(Y),'bo',lw=2)
xlim([-4,4])
ylabel(r"$|Y|$",size=16)
string = r"$\sin\left(1.25t\right)\times w(t)$"
title("Spectrum of "+string)
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y$",size=16)
12
xlabel(r"$\omega$",size=16)
grid(True)
#savefig("fig10-7.png")
show()

# Function definitions

def spectrum(function, T, N, Title, xLim = 10, lineWidth = 2):
    """
    Finds DFT and plots Magnitude and Phase spectra plots
    Inputs:
    function: function to calculate DFT for
    N: number of samples to be taken
    Title: Title for the plots
    xLim: xlimits for the plots
    lineWidth: linewidth for the plots
    """
    t=linspace(-T/2, T/2, N+1);t=t[:-1]
    dt=t[1]-t[0];fmax=1/dt
    y=function(t)
    y[0]=0 # the sample corresponding to -tmax should be set zero
    y=fftshift(y) # make y start with y(t=0)
    Y=fftshift(fft(y))/float(N)
    w=linspace(-pi*fmax,pi*fmax,N+1);w=w[:-1]
    figure()
    subplot(2,1,1)
    plot(w,abs(Y),lw=lineWidth)
    xlim([-xLim, xLim])
    ylabel(r"$|Y|$",size=16)
    title("Spectrum of "+Title)
    grid(True)
    subplot(2,1,2)
    plot(w,angle(Y),'ro',lw=lineWidth)
    xlim([-xLim, xLim])
    ylabel(r"Phase of $Y$",size=16)
    xlabel(r"$\omega$",size=16)
    grid(True)
    show()
    return Y

def spectrum_hamming(function, T, N, Title, xLim = 10, lineWidth = 2):
    """
    Finds DFT after applying Hamming window and plots Magnitude and Phase spectra plots
    Inputs:
    function: function to calculate DFT for
    N: number of samples to be taken
    Title: Title for the plots
    xLim: xlimits for the plots
    lineWidth: linewidth for the plots
    """
    t=linspace(-T/2, T/2, N+1);t=t[:-1]
    dt=t[1]-t[0];fmax=1/dt
    n=arange(N)
    wnd=fftshift(0.54+0.46*cos(2*pi*n/N))
    y=function(t)
    y=y*wnd
    y[0]=0 # the sample corresponding to -tmax should be set zero
    y=fftshift(y) # make y start with y(t=0)
    Y=fftshift(fft(y))/float(N)
    w=linspace(-pi*fmax,pi*fmax,N+1);w=w[:-1]
    figure()
    subplot(2,1,1)
    plot(w,abs(Y),lw=lineWidth)
    xlim([-xLim, xLim])
    ylabel(r"$|Y|$",size=16)
    title(r"Spectrum of "+Title)
    grid(True)
    subplot(2,1,2)
    plot(w,angle(Y),'ro',lw=lineWidth)
    xlim([-xLim, xLim])
    ylabel(r"Phase of $Y$",size=16)
    xlabel(r"$\omega$",size=16)
    grid(True)
    show()
    return Y

# Q2
t=linspace(-4*pi,4*pi,257);t=t[:-1]
def fn2(t):
    w0 = 0.86
    return cos(w0*t)**3

_ = spectrum(fn2, 8*pi, 256, r"$cos^3(0.86t)$ without Hamming window",)
_ = spectrum_hamming(fn2, 8*pi, 256, r"$cos^3(0.86t)$ with Hamming window",)

# Q3
random.seed(0)
w0 = 0.5 + random.random() # 0.5-1.5
d = random.random()*(2*pi) - pi # -pi to pi

T = 2*pi
N = 128
t=linspace(-T/2, T/2, N+1);t=t[:-1]
def fn3(t):
    return cos(w0*t+d)

Y_3 = spectrum_hamming(fn3, T, N, r"$cos(\omega_0t+\delta)$", 5)
# Weighted average
fmax = N/T
w = linspace(-pi*fmax,pi*fmax,N+1);w=w[:-1]
ii = where(w>=0)
w0_est = sum(w[ii]*abs(Y_3[ii])**3)/sum(abs(Y_3[ii])**3)
i = abs(w - w0_est).argmin()     
d_est = angle(Y_3[i])
print("Without noise:")
print("Actual w0:",w0)
print("Calculated w0:", w0_est)
print("Actual delta:",d)
print("Calculated delta:",d_est)

# Q4
def fn4(t):
    return cos(w0*t+d)+0.1*rand(N)

Y_4 = spectrum_hamming(fn4, T, N, r"$cos(\omega_0t+\delta)$", 5)
# Weighted average
noisy_w0_est = sum(w[ii]*abs(Y_4[ii])**3)/sum(abs(Y_4[ii])**3)
i = abs(w - noisy_w0_est).argmin()     
noisy_d_est = angle(Y_4[i])
print("With noise:")
print("Actual w0:",w0)
print("Calculated w0:", noisy_w0_est)
print("Actual delta:",d)
print("Calculated delta:",noisy_d_est)

T = 2*pi
N = 1024
t=linspace(-T/2, T/2, N+1);t=t[:-1]
def fn5(t):
    return cos(16*t*(1.5+t/(2*pi)))

_ = spectrum_hamming(fn5, T, N, r"$cos(16(1.5+\frac{t}{2\pi})t)$", 100)

T = 2*pi
N = 1024
t=linspace(-T/2, T/2, N+1);t=t[:-1]
split_time = split(t, 16)
DFT_mag = zeros((16, 64))
DFT_phase = zeros((16, 64))
dt=t[1]-t[0]; fmax=1/dt
n=arange(64)
wnd=fftshift(0.54+0.46*cos(2*pi*n/64))
for col in range(16):
    t_ = split_time[col]
    y=fn5(t_)
    y=y*wnd
    y[0]=0 # the sample corresponding to -tmax should be set zero
    y=fftshift(y) # make y start with y(t=0)
    Y=fftshift(fft(y))/64.0
    DFT_mag[col] = abs(Y)
    DFT_phase[col] = angle(Y)

t = t[::64]
w=linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
T, W = meshgrid(t, w)

fig1 = figure()
ax1 = p3.Axes3D(fig1)
fig1.add_axes(ax1)
title('Magnitude surface plot')
ax1.set_xlabel('t')
ax1.set_ylabel('w')
ax1.set_zlabel('Spectrum')
surf1 = ax1.plot_surface(W, T, DFT_mag.T, rstride=1, cstride=1, cmap=cm.jet)
fig1.colorbar(surf1)
show()
fig2 = figure()
ax2 = p3.Axes3D(fig2)
fig2.add_axes(ax2)
title('Phase surface plot')
ax2.set_xlabel('t')
ax2.set_ylabel('w')
ax2.set_zlabel('Phase')
surf2 = ax2.plot_surface(W, T, DFT_phase.T, rstride=1, cstride=1, cmap=cm.magma)
fig2.colorbar(surf2)
show()

