# file name: Recursive_Legendre_Temperature.py
#
# This script plots the thermal profile inside a sphere with one hemisphere at T0 and
# the other hemisphere at -T0, using the recursive definition of Legendre Polynomials.
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
# Parameters of the sphere
#
a = 2.0 # radius of the ring in meters
T0 = 100 # temperature of one hemisphere in Celsius
#
# Number of Legendre Polynomials to calculate, and values of xi = cos(theta)
#
N = 10 # first n+1 Legendre Polynomials (n>1)
r = 100 # number of values for theta between 0 and 2*pi, and for rho between 0 and a
theta = np.linspace(0, 2*np.pi, r) # theta values between 0 and 2*pi
xi = np.cos(theta) # xi = cos(theta)
#
# A function that computes a multiplicative coefficient for the potential
#
def Coeff(n):
    """Compute the multiplicative coefficient for the potential"""
    #
    Mul = 0.0
    for k in range(0, (n//2)+1):
        Mul += ((-1)**k)*math.factorial(2*n-2*k)/(math.factorial(k)*math.factorial(n-k)*math.factorial(n-2*k+1))
    Mul = ((2*n+1.0)/2.0**(n+1))*Mul
    #
    return Mul
#
# Legendre Polynomials are stored in a data frame such that column n is P_n(x)
#
mydata = pd.DataFrame(0.0, index=range(r), columns=range(N))
#
# Store P_0 and P_1
#
for i in range(r):
    mydata.at[i,0] = 1.0 # P_0(x) = 1
#
for i in range(r):
    mydata.at[i,1] = xi[i] # P_1(x) = xi
#
# Calculate values from P_2 to P_n using the recursive relation
#
for n in range(1,N-1):
    print(f"Calculating P_{n+1}(x)")
    for i in range(r):
        P1 = mydata.at[i,n-1] # P_(n-1)(x)
        P2 = mydata.at[i,n] # P_(n)(x)
        P3 = ((2*n+1)*xi[i]*P2 - n*P1)/(n+1) # P_(n+1)(x)
        mydata.at[i,n+1] = P3 # store P_(n+1)(x)
#
# Define the thermic field on plane y-z, as a function of theta and rho
#
rho = np.linspace(0, a, r)
T = np.zeros((r,r)) # temperature initialization
for i in range(r):
    for j in range(r):
        for n in range(1,N):
            T[i,j] += Coeff(n)*T0*Coeff(n)*mydata.at[j,n]*(rho[i]/a)**n
#
# Plot the temperature as a surface in 3D on plane Y-Z
#
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
Theta, Rho = np.meshgrid(theta, rho)
Z = Rho * np.cos(Theta) # convert to cartesian coordinates
Y = Rho * np.sin(Theta) # convert to cartesian coordinates
surf = ax.plot_surface(Y,Z,T,cmap='viridis')
#
# Axes, labels, view angle
#
ax.set_xlabel('y')
ax.set_ylabel('z')
ax.set_zlabel('Temperature (°C)')
ax.view_init(elev=30, azim=-45) # azim is rotation around vertical axis
plt.savefig('Temperature.jpg', dpi=300)
plt.close() # Close the current figure
#
# Plot isothermal curves on plane Y-Z
#
plt.figure(figsize=(6, 6))
Theta, Rho = np.meshgrid(theta, rho)
Y = Rho * np.sin(Theta)
Z = Rho * np.cos(Theta)
levels = 20  # number of isothermic curves
cs = plt.contour(Y, Z, T,levels=levels,cmap='viridis')
#
# Labels and title
#
plt.clabel(cs, inline=True, fontsize=8)
plt.xlabel('y')
plt.ylabel('z')
plt.title('Isothermal curves')
plt.axis('equal')
plt.savefig('Sphere_isothermal_curves.jpg', dpi=300)
plt.close() # Close the current figure
