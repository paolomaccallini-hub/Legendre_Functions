# file name: Recursive_Legendre_Expansion_1.py
#
# Derive and plot the expansion in Legendre polynomials of the step function
#
import os # functions for interacting with the operating system
import pandas as pd # dataframes similar to R
import matplotlib.pyplot as plt # a plotting library
import math # mathematical functions
import numpy as np # arrays and matrices
#
# Change the drectory, to the folder of the script
#
script_dir = os.path.dirname(os.path.abspath(__file__))
print(f"script_dir: {script_dir}")
os.chdir(script_dir)
#
# Number of Legendre Polynomials to calculate, and number of values of x 
#
N = 8 # first n+1 Legendre Polynomials (n>1)
r = 100 # number of values for x between -1 and 1
#
# Legendre Polynomials are stored in a data frame such that column n is P_n(x)
#
mydata = pd.DataFrame(0.0, index=range(r), columns=range(N))
#
# Store P_0 and P_1, 1 and x respectively
#
for i in range(r):
    mydata.at[i,0] = 1.0    # P_0(x) = 1
mydata[1] = np.linspace(-1, 1, r) # P_1(x) = x
#
# Calculate values from P_2 to P_n using the recursive relation
#
x = mydata[1] # x values
for n in range(1,N-1):
    print(f"Calculating P_{n+1}(x)")
    for i in range(r):
        P1 = mydata.at[i,n-1] # P_(n-1)(x)
        P2 = mydata.at[i,n] # P_(n)(x)
        P3 = ((2*n+1)*x[i]*P2 - n*P1)/(n+1) # P_(n+1)(x)
        mydata.at[i,n+1] = P3 # store P_(n+1)(x)
#
# Function f(x) to be expanded
#
def get_f(x):
    """Function to be expanded in Legendre Polynomials"""
    myfunction = pd.DataFrame(0.0, index=range(r), columns=range(1))
    for i in range(r):
        if x[i] < 0.0:
            myfunction.at[i,0] = 0.0
        else:
            myfunction.at[i,0] = 1.0
    return myfunction
#
# Calculate generalized Fourier coefficients a_n
#
mycoeff = pd.DataFrame(0.0, index=range(N), columns=range(1))
for i in range(N):
    Pn = mydata[i] # P_n(x)
    integrand = pd.DataFrame(0.0, index=range(r), columns=range(1))
    for j in range(r):
        integrand.at[j,0] = Pn[j]*get_f(x).at[j,0] # f(x)*P_n(x)
    mycoeff.at[i,0] = (2*i+1)/2 * np.trapezoid(integrand[0], x) # integral using trapezoidal rule
    print(f"a_{i} = {mycoeff.at[i,0]}")
#
# Build the interpolated function f_N(x)
#
f_N = pd.DataFrame(0.0, index=range(r), columns=range(1))
for i in range(N):
    Pn = mydata[i] # P_n(x)
    An = mycoeff.at[i,0] # a_n
    for j in range(r):
        f_N.at[j,0] += An*Pn[j] # f_N(x) += a_n*P_n(x)
#
# Plot 
#
plt.figure(figsize=(8, 5))
plt.plot(x, f_N[0], label='$f_N(x)$', color='black', linewidth=2) # interpolation
plt.plot(x, get_f(x)[0], label='$f(x)$', color='red', linewidth=2) # original function
for n in range(0,N): 
    Pn = mydata[n]
    An = mycoeff.at[n,0]
    plt.plot(x, An*Pn, label=f'$A_{n}P_{n}(x)$', linewidth=1, linestyle='--') # components
plt.xlabel('abscissa (x)')
plt.ylabel('$P_n(x)$')
plt.title('Expnsion in series of Legendre Polynomials')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('Legendre_Expansion_1.jpg', dpi=300)
plt.close() # Close the current figure