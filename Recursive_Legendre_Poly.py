# file name: Recursive_Legendre_Poly.py
#
# This script derives the values of the first n (specified by the user)
# Legendre Polynomials and plot them
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
N = 6 # first n+1 Legendre Polynomials (n>1)
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
# Plot 
#
plt.figure(figsize=(8, 5))
for n in range(0,N):
    Pn = mydata[n]
    plt.plot(x, Pn, label=f'$P_{n}(x)$')
plt.xlabel('abscissa (x)')
plt.ylabel('$P_n(x)$')
plt.title('Legendre Polynomials $P_n(x)$')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('Legendre_Polynomials.jpg', dpi=300)
plt.close() # Close the current figure