# 📡 Quantum Calculations: Solving the Schrödinger Equation Numerically

## 🔬 Introduction

Welcome to this repository, where we explore quantum mechanics numerically! Specifically, we solve the **time-independent Schrödinger equation** for a particle in a potential well using matrix methods. This is a fundamental problem in quantum mechanics, often used to model electrons in a quantum dot.

We use **Matlab** to construct and solve tridiagonal matrices that approximate the Schrödinger equation. Our goal is to compare numerical solutions with analytical results and investigate the behavior of quantum states under different potentials.

---

## 📜 Problem Formulation

The **time-independent Schrödinger equation** in one dimension is:

$$
-\frac{1}{2} \frac{d^2 \psi}{dx^2} + V(x) \psi(x) = E \psi(x)
$$

where:
- $\( \psi(x) \)$ is the wave function,
- $\( V(x) \)$ is the potential energy,
- $\( E \)$ is the energy eigenvalue.

### 🔹 Case 1: **Particle in a Box**
The simplest case is when $\( V(x) = 0 \)$ inside a box of width $\( a \)$, and infinite outside. The analytical solutions are:

$$
\psi_n(x) = \sqrt{\frac{2}{a}} \sin\left(\frac{n\pi x}{a}\right), \quad E_n = \frac{n^2 \pi^2}{2 a^2}
$$

where $\( n = 1, 2, 3, ... \)$ indexes the quantum states.

### 🔹 Case 2: **Harmonic-like Potential**
For a quantum dot system, the potential inside the well is modeled as:

$$
V(x) = \frac{(f x)^2}{2}
$$

where $\( f \)$ is a scaling factor. The goal is to find the **numerical eigenvalues and eigenfunctions** and compare them to first-order perturbation theory.

---

## 🛠️ Numerical Method

We discretize the Schrödinger equation using a **finite difference method**. The equation transforms into a matrix eigenvalue problem:

$$
H \psi = E \psi
$$

where $\( H \)$ is a **tridiagonal Hamiltonian matrix** with:
- **Main diagonal**: $\( \frac{1}{\Delta^2} + V_k \)$
- **Off-diagonal terms**: $\( -\frac{1}{2\Delta^2} \)$

Here, $\( \Delta \)$ is the discretization step size.

### 🔹 Matlab Implementation

Each script in this repository corresponds to a specific problem:

- **`kvantinl2_uppg2b.m`**: Computes the wave function for different values of \( f \) and visualizes how it changes.
- **`kvantinl2_uppg3a.m`**: Uses `fminsearch` to find the optimal \( f \) by minimizing the error between numerical and first-order perturbation theory results.
- **`kvantinl2_uppg3b.m`**: Investigates how first-order perturbation theory breaks down for large \( f \).

---

## 📊 Results & Discussion

### 🔹 Accuracy of Numerical Eigenvalues
We compare **numerical and analytical eigenvalues** for a particle in a box. As the number of grid points \( N \) increases, numerical results converge to analytical ones:

| Level \( n \) | Analytical \( E_n \) | Numerical \( E_n \) (N=10) | Numerical \( E_n \) (N=50) |
|--------------|----------------------|----------------------|----------------------|
| 1            | 4.9348                | 4.9014               | 4.9332               |
| 2            | 19.7392               | 19.2083              | 19.7143              |
| 3            | 44.4132               | 41.7619              | 44.2870              |
| ...          | ...                    | ...                  | ...                  |

For **higher energy levels**, numerical errors increase unless \( N \) is sufficiently large.

### 🔹 Effect of Potential Strength \( f \)
We plot the **ground-state wavefunction** for different values of \( f \):

- **For small \( f \)**: The wave function resembles that of a particle in a box.
- **For large \( f \)**: The wave function narrows and the energy increases significantly.

This follows from perturbation theory:

$$
E_1 \approx \frac{\pi^2}{2} + \frac{h^2 f^2 ( \pi - 6 )}{24 m_e \pi^2} \left(\frac{m_e a^2}{h^2} \right)
$$

For large \( f \), perturbation theory **fails** to capture the rapid energy increase.

---

## 📈 Figures

### 🔹 Numerical vs Analytical Wavefunction
<img src="wavefunction_comparison.png" width="500">

### 🔹 Effect of \( f \) on Ground State
<img src="wavefunction_vs_f.png" width="500">

---

## 🔮 Conclusion

- **Numerical methods** provide a powerful way to approximate quantum systems, especially when analytical solutions are not feasible.
- **Higher discretization points $ N $** improve accuracy but increase computational cost.
- **Perturbation theory is only valid for weak potentials**—for strong potentials, numerical solutions deviate significantly.
- **Future work:** Extend the model to **higher dimensions** and test different potentials!

---

## 📂 Repository Structure

