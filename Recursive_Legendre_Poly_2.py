# file name: Recursive_Legendre_Poly_2.py
#
# This script derives the values of the first n (specified by the user)
# Legendre Functions of the second kind and plot them
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
# Number of Legendre Functions of the second kind to calculate, and number of values of x 
#
N = 5 # first n+1 Legendre Functions of the second kind (n>1)
r = 100 # number of values for x between -1 and 1
eps = 1.0e-3 # small number to avoid singularities at x = +/- 1
x = np.linspace(-1+eps,1-eps, r) # x values between -1 and 1
#
# Legendre Functions of the second kind are stored in a data frame such that column n is Q_n(x)
#
mydata = pd.DataFrame(0.0, index=range(r), columns=range(N))
#
# Store Q_0 and Q_1, 1 and x respectively
#
for i in range(r):
    mydata.at[i,0] = 0.5*math.log((1+x[i])/(1-x[i])) # Q_0(x)
    mydata.at[i,1] = x[i]*mydata.at[i,0]-1 # Q_1(x)
#
# Calculate values from Q_2 to Q_n using the recursive relation
#
for n in range(1,N):
    print(f"Calculating Q_{n+1}(x)")
    for i in range(r):
        Q1 = mydata.at[i,n-1] # Q_(n-1)(x)
        Q2 = mydata.at[i,n] # Q_(n)(x)
        Q3 = ((2*n+1)*x[i]*Q2 - n*Q1)/(n+1) # Q_(n+1)(x)
        mydata.at[i,n+1] = Q3 # store Q_(n+1)(x)
#
# Plot 
#
plt.figure(figsize=(8, 5))
for n in range(0,N):
    Pn = mydata[n]
    plt.plot(x, Pn, label=f'$Q_{n}(x)$')
plt.xlabel('abscissa (x)')
plt.ylabel('$Q_n(x)$')
plt.title('Legendre Functions of the second kind $Q_n(x)$')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('Legendre_Functions of the second kind.jpg', dpi=300)
plt.close() # Close the current figure