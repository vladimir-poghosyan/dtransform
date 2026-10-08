# dtransform

**dtransform** is a symbolic Python library that implements the **Multivariable Differential Transform Method (DTM)**. It computes the *DTM spectrum* — a multivariate Taylor expansion — of symbolic expressions and supports arithmetic operations between spectra. Built on top of [SymPy](https://www.sympy.org/), it is designed to assist in symbolic modeling, semi-analytical computation, and solving PDEs/ODEs symbolically.

---

## ⚠️ Disclaimer

This project is part of ongoing research and may evolve as the work progresses. Interfaces and functionality are subject to change without notice. Use with caution in production or critical systems.

## 📘 What is Differential Transform Method?

The **Differential Transform Method (DTM)** is a semi-analytical technique that computes the Taylor series expansion of functions around a given point. It has applications in solving differential, integral, and functional equations.

Given a real-analytic multivariable function:

$$
f(x_1, x_2, \dots, x_n)
$$

its DTM coefficient at multi-index $\mathbf{k} = \left(k_1, k_2, \dots, k_n\right)$, centered at $\mathbf{a} = (a_1, a_2, \dots, a_n)$, with scaling constants $H_i > 0$, is defined as:

$$
F(k_1, \dots, k_n) =
\frac{H_1^{k_1} \cdots H_n^{k_n}}{k_1! \cdots k_n!}
\cdot
\left. \frac{\partial^{|\mathbf{k}|} f}{\partial x_1^{k_1} \cdots \partial x_n^{k_n}} \right|_{\mathbf{x} = \mathbf{a}}
$$

The inverse transform reconstructs the function as:

$$
f(x_1, \dots, x_n) \approx \sum_{k_1=0}^{N} \cdots \sum_{k_n=0}^{N}
F(k_1, \dots, k_n) \prod_{i=1}^n \left( \frac{x_i - a_i}{H_i} \right)^{k_i}
$$

This is the multivariable Taylor expansion expressed in normalized coordinates.

---

## ➕ Supported Operations

Let $\( A, B \)$ be DTM spectra, and $\( \lambda \in \mathbb{R} \)$ and $\( \boldsymbol\beta = (\beta_1, \dots, \beta_n) \)$ a multi-index of differentiation orders. Supported operations are:

- **Addition**:
  $\left(A + B\right)[\mathbf{K}] = A[\mathbf{K}] + B[\mathbf{K}]$

- **Subtraction**:
  $\left(A - B\right)[\mathbf{K}] = A[\mathbf{K}] - B[\mathbf{K}]$

- **Scalar Multiplication**:
  $\left(\lambda \cdot A\right)[\mathbf{K}] = \lambda \cdot A[\mathbf{K}]$

- **Cauchy Product (Convolution)**:
  $\left(A \cdot B\right)[\mathbf{K}] = \sum_{\mathbf{i} + \mathbf{j} = \mathbf{K}} A[\mathbf{i}] \cdot B[\mathbf{j}]$

- **Division**:

$$
C[\mathbf{K}] = \frac{1}{B[\mathbf{0}]}
\left( A[\mathbf{K}] - \sum_{\substack{\mathbf{i} + \mathbf{j} = \mathbf{K} \\ \mathbf{j} \ne \mathbf{K}}} B[\mathbf{i}] \cdot C[\mathbf{j}] \right)
$$

- **Differentiation**: (orders $\boldsymbol\beta$)

$$
G[\mathbf{K}] = \left(\prod_{i=1}^{n} \frac{(K_i + \beta_i)!}{K_i! \, H_i^{\beta_i}}\right)F[\mathbf{K} + \boldsymbol{\beta}]
$$

> ⚠️ Operations require matching variables, expansion centers, scaling constants, and orders.

---

## 🧠 Features

- Symbolic DTM spectrum computation
- Multivariable expansions with arbitrary centers and scalings
- Symbolic inverse transformation
- Supports `+`, `-`, `*`, `/`, and scalar operations
- Spectrum-domain differentiation
- Works seamlessly with SymPy
- Easily extensible

---

## 📦 Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/vladimir-poghosyan/dtransform.git
cd dtransform
python3 -m pip install .
```

## 🛠️ API

```Python
Spectrum(
    expr: str,
    order: int = 4,
    center: dict[str, int | float] = None,
    scaling: dict[str, int | float] = None,
    **kwargs: int | float
)
```

Parameters

    expr:
        A string representing the symbolic expression (e.g., "x + y").

    order:
        Maximum order of expansion (default: 4).

    center (optional):
        Dictionary specifying the expansion center for each variable, e.g., {"x": 1, "y": 2}.
        If omitted, all variables default to 0.

    scaling (optional):
        Dictionary specifying normalization factors (scaling constants) for each variable.
        If omitted, all variables default to a scaling of 1.

    **kwargs:
        Optional variable center values supplied directly as keyword arguments.
        These override corresponding values in center.

Methods

    inverse() — Reconstructs the symbolic function from its DTM spectrum.

    diff(variable: str | None = None, /, **orders: int) — Computes the DTM spectrum of a derivative.

    display_coefficients() — Displays non-zero coefficients of the transformation spectrum.

    Supports arithmetic operators: +, -, *, /, including scalar multiplication and division.

## 🚀 Usage Example

```Python
from dtransform import Spectrum


# Define symbolic expressions and expansion center
expr1 = "x + y"
expr2 = "1 + x * y"
center = {'x': 1, 'y': 2}

# Create DTM spectra
s1 = Spectrum(expr1, order=4, center=center)
s2 = Spectrum(expr2, order=4, center=center)

# Perform operations
sum_spectrum = s1 + s2
prod_spectrum = s1 * s2
quot_spectrum = s1 / s2

# Differentiate spectra
dx_spectrum = s1.diff("x")
dxy_spectrum = s2.diff(x=1, y=1)

# Reconstruct and print inverse transforms
print("Reconstructed Sum:      ", sum_spectrum.inverse())
print("Reconstructed Product:  ", prod_spectrum.inverse())
print("Reconstructed Quotient: ", quot_spectrum.inverse())
print("Reconstructed d/dx: ", dx_spectrum.inverse())
print("Reconstructed d²/dxdy: ", dxy_spectrum.inverse())
```

Output (simplified):

```
Reconstructed Sum:      x + y + x*y + 1
Reconstructed Product:  x**2*y + x*y**2 + x + y
Reconstructed Quotient: 4*x**3*y**3/729 - 5*x**3*y**2/243 - 32*x**3*y/243 + 112*x**3/729 - x**2*y**3/81 + 7*x**2*y/9 - 64*x**2/81 - 4*x*y**3/243 + 23*x*y**2/81 - 160*x*y/81 + 419*x/243 + 17*y**3/729 - 64*y**2/243 + 323*y/243 - 64/729
Reconstructed d/dx: 1
Reconstructed d²/dxdy: 1
```

## 📚 References

1. **Pukhov, A. A.**, Differential Transformation Method and Its Application in Mechanics and Control Theory. Moscow State University, 1971. (in Russian)
2. **Zhou, J. K.**, Differential Transformation and Its Applications for Electrical Circuits. Huazhong University Press, 1986.
3. **Simonyan S. H., Avetisyan A. G.**, Applied theory of differential transforms. Chartaraget, Yerevan, 2010.
4. **Poghosyan V. A., Avetisyan A. G.**, The multidimensional differential analog of Berkowitz’s algorithm. Proceedings of the RA National Academy of Science and NPUA, Yerevan, 2017.
5. **Bervillier C.**, Status of the differential transformation method. Université François Rabelais, 2012.
6. **Abdel-Halim Hassan**, Application to differential transformation method for solving systems of differential equations. Zagazig University, 2007.
7. **Arikoglu, A. and Ozkol, I.**, Solution of differential-difference equations by using differential transform method, Applied Mathematics and Computation, 2005.
8. **Hassan, I. H.**, Comparison between differential transformation method and Adomian decomposition method, Applied Mathematics and Computation, 2008.

## 📄 License

This project is licensed under the MIT License. See LICENSE for details.
