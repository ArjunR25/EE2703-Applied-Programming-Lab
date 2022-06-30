'''
Laplace Equation
- Solves the resistor problem using Laplace equation for potential
- Finds the potential matrix subject to user inputs Nx, Ny, radius and Niter
- Obtains the current distribution in the resistor using this potential matrix

Mar 08, 2022

By Arjun R
EE20B016
'''

from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
from scipy.linalg import lstsq
from sys import argv, exit

'''Taking input'''
if len(argv) == 5:
    Nx = int(argv[1]) # size along x
    Ny = int(argv[2]) # size along y
    radius = int(argv[3]) # radius of central lead
    Niter = int(argv[4]) # number of iterations to perform
    if Nx != Ny:
        print("Error! Nx doesn't match Ny")
        exit()
    elif radius > Nx:
        print("Error! Entered radius should be less than Nx and Ny values")
        exit()

    if Niter <= 500:
        print("Error! Entered Niter value should be greater than 500")
        
else:
    #default
    Nx = 25
    Ny = 25
    radius = 8
    Niter = 1500

try:
    '''Function definitions'''
    def update_phi(phi): # Function to update the potential matrix
        phi[1:-1,1:-1] = 0.25*(phi[1:-1,0:-2]+phi[1:-1,2:]+phi[0:-2,1:-1]+phi[2:,1:-1])
        phi[:,0] = phi[:,1] # left boundary 
        phi[:,-1] = phi[:,-2] # right boundary 
        phi[0,:] = phi[1,:] # top boundary
        # bottom boundary is not updated at all, hence stays zero

    def exponential(A, B): # Calculates the exponential function given A and B parameters
        x = arange(Niter)
        return A*exp(B*x)

    def get_cum_error(A, B, N): # Returns the cumulative error acc. to formula
        return -(A/B)*exp(B*(N+0.5))

    '''Defining the potential matrix'''
    phi = np.zeros((Ny, Nx)) # Initializing the potential matrix

    x = np.linspace(-0.5, 0.5, Nx)
    y = np.linspace(-0.5, 0.5, Ny)
    X, Y = meshgrid(x, y)

    r = radius/(Nx-1) # e.g., 0.5*(8/12) = 0.35

    ii = where(X*X + Y*Y < (r+0.01)**2)  
    x_lead, y_lead = ii # returns the x and y-coordinates covered by the central lead

    phi[ii] = 1.0 # Central lead is held at 1 V

    # First Contour Plot
    figure(1)
    cp = contourf(X, Y[::-1], phi, cmap = cm.plasma, levels = 100)
    cbar = colorbar(cp)
    plot((x_lead-(Nx-1)/2)/(Nx-1), (y_lead-(Ny-1)/2)/(Ny-1), 'ro', label = 'lead points')
    xlabel('X')
    ylabel('Y')
    legend()
    title('Contour plot for potential')
    show()

    errors = zeros(Niter) # For storing max. error in each iteration
    for k in range(Niter):
        oldphi=phi.copy()
        update_phi(phi)
        phi[ii] = 1.0 # over-written values at the central lead

        errors[k] = (abs(phi-oldphi)).max()

    '''Finding fit'''
    M = np.zeros((Niter, 2)) # 
    M[:, 0] = 1
    p = np.zeros(Niter)
    for i in range(Niter):
        M[i, 1] = i
        p[i] = log(errors[i])

    x1 = lstsq(M, p)[0] # for whole list of errors
    x2 = lstsq(M[500:,:], p[500:])[0] # after 500 iterations
    # returns logA, B
    # rcond = None

    #print(x1, x2)

    figure(2)
    semilogy(errors, label='errors', color = 'blue')
    semilogy(arange(0, Niter, 50), errors[::50], 'ro')
    semilogy(exponential(exp(x1[0]), x1[1]), label='fit 1', color = 'green')
    semilogy(exponential(exp(x2[0]), x2[1]), label='fit 2', color = 'black')
    xlabel('iterations')
    ylabel('Error')
    legend()
    grid()
    title('Best fit plot for error (semilogy scale)')
    show()

    figure(3)
    loglog(errors, label='error', color = 'blue')
    loglog(arange(0, Niter, 50), errors[::50], 'ro')
    loglog(exponential(exp(x1[0]), x1[1]), label='fit 1', color = 'green')
    loglog(exponential(exp(x2[0]), x2[1]), label='fit 2', color = 'black')
    xlabel('iterations')
    ylabel('Error')
    legend()
    grid()
    title('Best fit plot for error (log-log scale)')
    show()

    '''Cumulative Error'''
    cum_error = get_cum_error(exp(x2[0]), x2[1], np.arange(Niter))

    figure(4)
    loglog(np.arange(100, Niter, 100), cum_error[100::100], 'ro')
    xlabel('iterations')
    ylabel('Maximum cumulative error')
    grid()
    title('Log-log plot for cumulative error')
    show()

    '''3D Surface plot'''
    fig5 = figure(5) # open a new figure
    ax = p3.Axes3D(fig5, auto_add_to_figure=False) # Axes3D is the means to do a surface plot
    fig5.add_axes(ax)
    title('The 3-D surface plot of the potential')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Potential')
    surf = ax.plot_surface(X, Y[::-1], phi, rstride=1, cstride=1, cmap=cm.jet) #phi.T
    show()

    '''Final contour plot'''
    figure(6)
    cp = contourf(X, Y[::-1], phi, cmap = cm.plasma, levels = 100)
    cbar = colorbar(cp)
    plot((x_lead-(Nx-1)/2)/(Nx-1), (y_lead-(Ny-1)/2)/(Ny-1), 'ro', label = 'lead points')
    title('The contour plot of the potential')
    xlabel('X')
    ylabel('Y')
    legend()
    show()

    '''Finding the current distribution'''
    Jx = 0.5*(phi[1:-1,0:-2]-phi[1:-1,2:]) # left - right
    Jy = 0.5*(phi[0:-2,1:-1]-phi[2:,1:-1]) # top - bottom

    figure(7)
    quiver(X[1:-1,1:-1], Y[-1:1:-1,-1:1:-1], Jx, -Jy, label = r'$\vec{J}$')
    #quiver(X[1:-1,1:-1], Y[1:-1,1:-1], Jx[::-1,:], Jy[::-1,:])
    plot((x_lead-(Nx-1)/2)/(Nx-1), (y_lead-(Ny-1)/2)/(Ny-1), 'ro', label = 'Lead points')
    xlabel('X')
    ylabel('Y')
    legend()
    title('The vector plot of the current flow')
    axis('square')
    show()
except:
    print("Error!")
    exit()
