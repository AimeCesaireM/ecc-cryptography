# Elliptic Curve Cryptography

## Description

This repository contains the implementation of elliptic curve cryptography algorithms in Python, including several cryptographic functions and methods for working with elliptic curves over finite fields. The code is part of a project for the course "Cryptographie à base de courbes elliptiques" at Sorbonne Université, LU3IN024 (2021-2022).

The project includes key elliptic curve cryptographic operations, such as point addition, scalar multiplication, and factorization, along with additional utility functions for curve analysis and cryptographic computations.

## Key Features

- **Elliptic Curve Validation**: Functions to verify whether a curve is elliptic and to check if a point belongs to the curve.
- **Point Operations**: Implementations of elliptic curve point addition, point negation, and scalar multiplication.
- **Curve Analysis**: Functions to compute the cardinality of elliptic curves, generate a list of points on a curve, and verify curve properties.
- **Mathematical Utilities**: Includes modular exponentiation, modular inverse, Legendre symbol calculation, and more.

## Prerequisites

To run the code, you need to have Python 3 installed along with the following packages:

- `matplotlib` for plotting graphs of curve cardinalities.

You can install the required packages using `pip`:

```bash
pip install matplotlib
```

## Functions

### Utility Functions

- `exp(a, N, p)`: Performs fast modular exponentiation, calculating \(a^N \mod p\).
- `factor(n)`: Returns the prime factorization of a number \(n\).
- `inv_mod(x, p)`: Returns the modular inverse of \(x\) modulo \(p\).
- `racine_carree(a, p)`: Computes the square root of \(a\) modulo \(p\) when \(p \equiv 3 \mod 4\).

### Elliptic Curve Operations

- `est_elliptique(E)`: Checks if a curve \(E = (p, a, b)\) is elliptic.
- `point_sur_courbe(P, E)`: Checks if a point \(P = (x, y)\) is on the elliptic curve \(E\).
- `addition(P1, P2, E)`: Performs elliptic curve point addition on points \(P1\) and \(P2\) on curve \(E\).
- `multiplication_scalaire(k, P, E)`: Performs scalar multiplication of a point \(P\) on the elliptic curve \(E\).

### Curve Properties

- `cardinal(E)`: Returns the cardinality of the group of points on the elliptic curve \(E\).
- `liste_points(E)`: Returns a list of points on the elliptic curve \(E\).
- `verif_liste_points(E)`: Verifies that all points in `liste_points(E)` satisfy the curve equation.

### Random Point Generation

- `point_aleatoire_naif(E)`: Generates a random point on the elliptic curve \(E\) using a naive approach.
- `point_aleatoire(E)`: Generates a random point on the elliptic curve \(E\) using a more efficient method.

### Curve Visualization

- `dessine_graphe(p)`: Plots a histogram of the distribution of the cardinalities of elliptic curves over the finite field \(F_p\).

## Example Usage
### Example 1
```python
#
p = 248301763022729027652019747568375012323
E = (p, 1, 0)
N = 248301763022729027652019747568375012324
facteurs = [(2, 2), (62075440755682256913004936892093753081, 1)]

pointDHM = diffie_hellman(E, N, facteurs)

```
### Example 2
```python
# Define the curve parameters
p = 360040014289779780338359
a = 117235701958358085919867
b = 18575864837248358617992
E = (p, a, b)

# Check if the curve is elliptic
if est_elliptique(E):
    print("The curve is elliptic.")

# Generate a random point on the curve
point = point_aleatoire(E)
print(f"Random point on the curve: {point}")

# Calculate the cardinality of the elliptic curve
cardinality = cardinal(E)
print(f"Cardinality of the curve: {cardinality}")
```

## License

This project is licensed under the MIT License
