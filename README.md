# Heat transfer equation FEM solver
Differential Equations course at AGH UST - 2024/2025.

## Problem

Find the function $u \colon [0, 2] \to \mathbb{R}$  satisfying the following equation:

$$
\begin{align*} 
-\frac{d}{dx} \left( k(x) \frac{du(x)}{dx} \right) = 100x^2, \quad x \in [0, 2], \\
\\
\\
u(2) = -20, \\
\\
\\
\frac{du(0)}{dx} + u(0) = 20, \\
\\
\\
k(x) =
\begin{cases}
1 & \text{for } x \in [0, 1], \\
2x & \text{for } x \in (1, 2].
\end{cases}
\end{align*}
$$

Weak formulation for this problem is derived in [`FEM-weak-formulation.pdf`](FEM-weak-formulation.pdf).

For a maximum grade, integration must be performed numerically ([**Gauss–Legendre quadrature**](https://en.wikipedia.org/wiki/Gauss–Legendre_quadrature) method) and the number of elements must be defined by the user.

Grade - 50/50pts.
