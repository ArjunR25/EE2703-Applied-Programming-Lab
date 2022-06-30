'''
CURRENT DISTRIBUTION IN A HALFWAVE DIPOLE ANTENNA
- Takes in dipole parameters as input
- Finds the exact current distribution
- Compares it with the assumed current distribution

Mar 27, 2022

By Arjun R
EE20B016
'''
from pylab import *

### 1 ###
N = 100 # number of sections in each half-section
l = 0.5 # quarter wavelength
Im = 1.0 # current injected at feed point
a = 0.01 # radius
dz = l/N # spacing b/w current samples
mu0 = 4e-7 * pi # permeability of free space
k = pi/(2*l) # wavenumber

z = arange(-N, N+1)*dz # ranges from -N*dz to N*dz
u = delete(z, [0, N, 2*N]) # removes locations where current is known
#I = zeros(2*N+1)
#J = zeros(2*N-2)

### 2 ###
def return_M(N, a): 
    """Returns M matrix where H = M*J
    Inputs:
    N = number of sections
    a = radius at observing point"""
    M = identity(2*N-2)/(2*pi*a) # identity(N) returns identity matrix of size N
    return M

### 3 ###
zj, zi = meshgrid(z, z) # creates a rectangular grid out of x and y arrays
uj, ui = meshgrid(u, u) 
Rz = array(sqrt(a**2 + (zi-zj)**2)) # distance matrix
Ru = array(sqrt(a**2 + (ui-uj)**2)) # matrix containing distances to unknown currents
RiN = delete(Rz[:, N], [0, N, 2*N]) # size=2N-1; matrix with distance to center element Im 

# Finding P matrices
P = mu0*(1/Ru)*exp(-(1j)*k*Ru)*dz/(4*pi) 
PB = mu0*(1/RiN)*exp(-(1j)*k*RiN)*dz/(4*pi) 

### 4 ###
def vec_complex(a, b): 
    """
    returns complex number given real and imaginary values"""
    return complex(a, b)

vec_complex = vectorize(vec_complex) # vectorizes function for taking in arrays

# Finding Q matrices
Q = P*a*vec_complex(1/(Ru**2), k*(1/Ru))/mu0
QB = PB*a*vec_complex(1/(RiN**2), k*(1/RiN))/mu0

### 5 ###
M = return_M(N, a) 
J = Im * dot(inv(M-Q), QB) # Calculates J array
I = J.copy()
I = insert(I, 0, 0)
I = insert(I, N, Im)
I = insert(I, 2*N, 0)
# calculates final I array after inserting center and boundary values 

I_assumed = Im*sin(k*(l+abs(z))) 
# calculates initially assumed sinusoidally distributed current

# Plotting I vs z
fig = figure(0)
plot(z, abs(I), 'b', label = 'Calculated value')
plot(z, I_assumed, 'r', label = 'Assumed value')
xlabel('z')
ylabel('Current')
title('Current vs z')
legend()
grid()
show()

print("M:\n",M.round(2))
print("z:\n",z.round(2))
print("u:\n",u.round(2))
print("Rz:\n",Rz.round(2))
print("Ru:\n",Ru.round(2))
print("RiN:\n",RiN.round(2))
print("P*((10)*8):\n",(P*1e8).round(2))
print("PB*((10)*8): \n",(PB*1e8).round(2))
print("Q*((10)*5):\n",(Q*1e5).round(2))
print("QB*((10)*5):\n",(QB*1e5).round(2))
print("I:\n",I.round(2))
print("J * 10*3 :\n",(J*1e3).round(2))
