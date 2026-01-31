# Repository documentation

**Full document (PDF):**  
[View the PDF](Legendre%20Functions.pdf)

**Example: recursion formula**

The first N Legendre polynomials can be derived from $P_0$ and $P_1$ using the recursion formula:

$$
P_{n+1}(x)=\frac{2n+1}{n+1}xP_{n}(x)-\frac{n}{n+1}P_{n-1}(x)
$$

The script `Recursive_Legendre_Expansion_1.py` calculates the first N Legendre polynomials using the recursion formula and plots them, as in the example below.

![Legendre_Polynomials](https://github.com/user-attachments/assets/20daf15d-2e38-449d-9d2f-433757d9b67e)

A formally identical recursion formula holds for Legendre functions of the second type. The script `Recursive_Legendre_Expansion_2.py` calculates the first N Legendre functions of the second type using the recursion formula and plots them, as in the example below.

![Legendre_Functions of the second kind](https://github.com/user-attachments/assets/f74f09d7-d0da-4912-9027-0ed2d67883c7)

**Example: expansion in Legendre polynomials**

Any function generally continuous in [0,1] can be expanded in a series of Legendre polynomials. `Recursive_Legendre_Expansion_2.py` calculates the expansion for the step function (first N terms) and plots the approximated function and the N components, as shown below.

![Legendre_Expansion_2](https://github.com/user-attachments/assets/8a25e21c-0d6d-4193-a669-73f87cc9f5ec)

**Example: gravitational potential of a ring**

Legendre polynomials $P_n$ and functions of the second kind $Q_n$ allow us to express the general integral of the Laplace equation in spherical coordinates, when the problem is independent of the azimuth:

$$
\nabla^2 f =
\frac{1}{r^2}\frac{\partial}{\partial r}
\left(r^2 \frac{\partial f}{\partial r}\right)
+
\frac{1}{r^2 \sin\theta}\frac{\partial}{\partial \theta}
\left(\sin\theta \frac{\partial f}{\partial \theta}\right)
= 0
$$

The solution is of the following type:

$$
f(\rho,\theta)=\left[C_1 P_n(\theta) + C_2 Q_n(\theta)\right]\left[C_3 \rho^n + C_4 \rho^{-(n+1)}\right]
$$

with $C_1, C_2, C_3, C_4$ arbitrary constants. One application of this result is the calculation of the gravitational potential of a ring in the space around it. `Recursive_Legendre_Gravity.py` finds the solution for the potential of a ring placed on the plane xy, witha given radius _a_ and a given linear density $\sigma$ using the first N Legendre polynomials (the Legendre functions of the second type have no physical meaning in this case, since they are divergent) and plots it as a surface on the plane yz, along with the euipotential lines, as shown below.

![Ring_potential_equipotential](https://github.com/user-attachments/assets/3f97c17d-069c-4399-b177-627bdb471625)
