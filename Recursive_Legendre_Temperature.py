# file name: Recursive_Legendre_Gravity.py
#
# This script plots the graviational potential of a ring mass
#
import os # functions for interacting with the operating system
import pandas as pd # dataframes similar to R
import matplotlib.pyplot as plt # a plotting library
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (needed for 3D)
import math # mathematical functions
import numpy as np # arrays and matrices
#
# Change the drectory, to the folder of the script
#
script_dir = os.path.dirname(os.path.abspath(__file__))
print(f"script_dir: {script_dir}")
os.chdir(script_dir)
#
# Parameters of the ring
#
a = 2.0 # radius of the ring in meters
sigma = 1.0 # mass per unit length of the ring in kg/m
G = 6.67430e-11 # gravitational constant in m^3 kg^-1 s^-2
#
# Number of Legendre Polynomials to calculate, and values of xi = cos(theta)
#
N = 50 # first n+1 Legendre Polynomials (n>1)
r = 100 # number of values for theta between 0 and 2*pi, and for rho between 0 and a
theta = np.linspace(0, math.pi, r) # theta values between 0 and 2*pi
xi = np.cos(theta) # xi = cos(theta)
#
# A function that computes a multiplicative coefficient for the potential
#
def Coeff(n):
    """Compute the multiplicative coefficient for the potential"""
    Mul = ((-1)**n)*math.factorial(2*n)/(2**(2*n)*(math.factorial(n))**2)
    return Mul
#
# Legendre Polynomials are stored in a data frame such that column n is P_n(x)
#
mydata = pd.DataFrame(0.0, index=range(r), columns=range(2*N))
#
# Store P_0 and P_1
#
for i in range(r):
    mydata.at[i,0] = 1.0    # P_0(x) = 1
#
for i in range(r):
    mydata.at[i,1] = xi[i]          # P_1(x) = xi
#
# Calculate values from P_2 to P_n using the recursive relation
#
for n in range(1,2*N-1):
    print(f"Calculating P_{n+1}(x)")
    for i in range(r):
        P1 = mydata.at[i,n-1] # P_(n-1)(x)
        P2 = mydata.at[i,n] # P_(n)(x)
        P3 = ((2*n+1)*xi[i]*P2 - n*P1)/(n+1) # P_(n+1)(x)
        mydata.at[i,n+1] = P3 # store P_(n+1)(x)
#
# Define the gravitational potential of a ring on plane y-z, as a function of theta and rho
#
# Case rho <= a
#
rho1 = np.linspace(0, a, r)
V1 = np.zeros((r,r)) # potential initialization
for i in range(r):
    for j in range(r):
        for n in range(N):
            V1[i,j] += 2.0*math.pi*G*sigma*Coeff(n)*mydata.at[j,2*n]*(rho1[i]/a)**(2*n)
#
# Case rho > a
#
rho2 = np.linspace(a, 5*a, r)
V2 = np.zeros((r,r)) # potential initialization
for i in range(r):
    for j in range(r):
        for n in range(N):
            V2[i,j] += 2.0*math.pi*G*sigma*Coeff(n)*mydata.at[j,2*n]*(a/rho2[i])**(2*n+1)
#
# Plot the potential as a surface in 3D on plane Y-Z
#
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#
# Set mask for singularity at rho=a, theta=pi/2
#
delta_Rho = a/50
delta_Theta = math.pi/12
#
# Plot V1
#
Theta, Rho = np.meshgrid(theta, rho1) # grid for theta (column) and rho (row)
# mask = ((np.abs(Rho-a)<delta_Rho)&(np.abs(Theta-math.pi/2)<delta_Theta)) # mask for singularity at rho=a, theta=pi/2
V1[mask] = np.nan # set singularity to NaN
Z = Rho * np.cos(Theta) # convert to cartesian coordinates
Y = Rho * np.sin(Theta) # convert to cartesian coordinates
surf = ax.plot_surface(Y,Z,V1, cmap='viridis')
#
# Plot V2
#
Theta, Rho = np.meshgrid(theta, rho2) # grid for theta (column) and rho (row)
# mask = ((np.abs(Rho-a)<delta_Rho)&(np.abs(Theta-math.pi/2)<delta_Theta)) # mask for singularity at rho=a, theta=pi/2
V2[mask] = np.nan # set singularity to NaN
Z = Rho * np.cos(Theta) # convert to cartesian coordinates
Y = Rho * np.sin(Theta) # convert to cartesian coordinates
surf = ax.plot_surface(Y,Z,V2, cmap='viridis')
#
# Axes, labels, view angle
#
ax.set_xlabel('y')
ax.set_ylabel('z')
ax.set_zlabel('V')
ax.view_init(elev=30, azim=30) # azim is rotation around vertical axis
plt.savefig('Ring_potential.jpg', dpi=300)
plt.close() # Close the current figure
#
# Plot equipotential curves on plane Y-Z
#
plt.figure(figsize=(6, 6))
#
# V1
#
Theta, Rho = np.meshgrid(theta, rho1)
Y = Rho * np.sin(Theta)
Z = Rho * np.cos(Theta)
levels = 20  # number of equipotential curves
cs = plt.contour(Y, Z, V1,levels=levels,cmap='viridis')
#
# V2
#
Theta, Rho = np.meshgrid(theta, rho2)
Y = Rho * np.sin(Theta)
Z = Rho * np.cos(Theta)
levels = 20  # number of equipotential curves
cs = plt.contour(Y, Z, V2,levels=levels,cmap='viridis')
#
# Labels and title
#
plt.clabel(cs, inline=True, fontsize=8)
plt.xlabel('y')
plt.ylabel('z')
plt.title('equipotential curves')
plt.axis('equal')
plt.savefig('Ring_equipotential_curves.jpg', dpi=300)
plt.close() # Close the current figure