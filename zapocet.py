import configparser
import matplotlib.pyplot as plt
import numpy as np


# Set values from config file
config = configparser.ConfigParser()
config.read('Zapocet\\advection.ini')

# Store values into variables
nx = int(config.get('DEFAULT', 'nx'))
nt = int(config.get('DEFAULT', 'nt'))
dx = int(config.get('DEFAULT', 'dx'))
dt = int(config.get('DEFAULT', 'dt'))
c = int(config.get('DEFAULT', 'c'))
plot_every = int(config.get('DEFAULT', 'plot_every'))

# Load initial conditions for phi in t=0
ic_loaded = np.genfromtxt('Zapocet\\ic.dat')
ic = np.full(len(ic_loaded)+2, 0.0)
# Create an array for initial conditions with the boundary condition in the assignment
ic[0] = ic_loaded[-1]
ic[-1] = ic_loaded[0]
ic[1:-1] = ic_loaded


# x-grid points as an array
x_array_loaded = np.arange(nx)
x_array = np.full(len(x_array_loaded)+2, 0.0)
x_array[0] = -1
x_array[-1] = x_array_loaded[-1] + 1
x_array[1:-1] = x_array_loaded

# Time steps as an array
t_array = np. arange(nt)

# Creates a matrix with values of phi where going down the rows means moving in time and going left to right in columns means moving in x-direction
phi_matrix = np.full((len(t_array), len(ic)), 0.0)

# Boundary conditions
phi_matrix[0] = ic
phi_matrix[-1] = ic


# Using 'frog-leap', values for phi in time t+1 are calculated
for i in range(len(phi_matrix)-1): # this is time
    for j in range(1, len(phi_matrix[i])-1): # this is x coordinates
        phi_matrix[i+1][j] = phi_matrix[i-1][j] - (c * dt * ((phi_matrix[i][j+1] - phi_matrix[i][j-1])/dx))
    phi_matrix[i+1][0] = phi_matrix[i+1][-2]
    phi_matrix[i+1][-1] = phi_matrix[i+1][1]

# Matplotlib stuff for plotting a series of graphs
# Creates 100 frames, which could be then used for making a gif
frames = 100
for i in range(frames):
    fig, ax = plt.subplots()
    ax.set_title("t = " + str(i*dt) + "s", loc = 'left')
    ax.set_xlabel("x")
    ax.set_ylabel("Ψ")
    ax.plot(x_array[1:-1]*dx, phi_matrix[i*plot_every][1:-1], color='#c736c4', label='Advanced numerical solution')
    ax.legend(loc = 'upper right')
    plt.ylim(-2, 2)
    plt.xlim(x_array[1]*dx, x_array[-1]*dx)
    ax.grid(linestyle='-', linewidth=1)
    plt.savefig("fig" + str(i))

# Code by Rudolf Pardubicky
