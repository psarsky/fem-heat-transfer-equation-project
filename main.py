from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt

def k(x: float) -> float:
    return 1 if x <= 1 else 2 * x

def element_integrals(x_i: float, x_ip1: float, h: float) -> Tuple[np.ndarray, np.ndarray]:
    B_local: np.ndarray = np.zeros((2, 2))
    L_local: np.ndarray = np.zeros(2)
    quadrature_points: np.ndarray = np.array([-1, 1]) / np.sqrt(3)
    midpoint: float = np.mean([x_i, x_ip1])
    de_dx: np.ndarray = np.array([-1, 1]) / h

    for xi in quadrature_points:
        x: float = midpoint + h * xi / 2
        e: np.ndarray = np.array([1 - xi, 1 + xi]) / 2
        for i in range(2):
            L_local[i] += 100 * x**2 * e[i] * h / 2
            for j in range(2):
                B_local[i, j] += k(x) * de_dx[i] * de_dx[j] * h / 2

    return B_local, L_local

def assemble_global_matrices(N: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    nodes: np.ndarray = np.linspace(0, 2, N + 1)
    B_global: np.ndarray = np.zeros((N + 1, N + 1))
    L_global: np.ndarray = np.zeros(N + 1)

    for element in range(N):
        B_local, L_local = element_integrals(nodes[element], nodes[element + 1], 2 / N)
        for i in range(2):
            L_global[element + i] += L_local[i]
            for j in range(2):
                B_global[element + i, element + j] += B_local[i, j]

    B_global[0, 0] += 1
    L_global[0] -= 20

    B_global[-1, :] = 0
    B_global[-1, -1] = 1
    L_global[-1] = -20

    return B_global, L_global, nodes

def solve_and_plot(n: int) -> None:
    B_global, L_global, nodes = assemble_global_matrices(n)
    u: np.ndarray = np.linalg.solve(B_global, L_global)

    print(B_global)
    print(L_global)

    plt.plot(nodes, u, "o--", label=f"FEM Solution (N={n} elements)")
    plt.title("Heat Transfer Equation FEM Solution")
    plt.xlabel("Position x")
    plt.ylabel("Temperature u(x)")
    plt.grid()
    plt.legend()
    plt.show()

if __name__ == "__main__":
    elements: int = int(input("Provide number of elements: "))
    solve_and_plot(elements)
