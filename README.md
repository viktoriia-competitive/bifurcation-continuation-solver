# Bifurcation Continuation Solver

![Python](https://img.shields.io/badge/Python-3.x-blue)
![Dependencies](https://img.shields.io/badge/Dependencies-Standard%20Library-green)
![Topic](https://img.shields.io/badge/Topic-Dynamical%20Systems-purple)
![Numerics](https://img.shields.io/badge/Numerics-Newton%20%2B%20Continuation-orange)

A compact numerical continuation script for tracking fixed points of parameter-dependent quadratic maps and detecting when an eigenvalue crosses a selected target value.

## Overview

This project studies fixed points of a map depending on a continuation parameter.  
The script treats the unknown vector as:

```text
p = (x_1, x_2, ..., x_n, parameter)
```

and solves the fixed-point equations:

```text
F_i(p) = f_i(p) - x_i = 0,   i = 1, ..., n
```

where each `f_i` is represented as a quadratic polynomial in all `n + 1` variables.

The solution set is generically a one-dimensional curve in `(state, parameter)` space. The script follows this curve using a predictor-corrector method and detects when an eigenvalue of the Jacobian crosses a target value.

## What the script does

1. Reads a fixed-point continuation problem from standard input.
2. Infers the dimension `n` from the amount of polynomial coefficient data.
3. Builds the quadratic monomial basis.
4. Evaluates the map and its Jacobian.
5. Computes a tangent direction to the fixed-point curve.
6. Advances along the curve using a small predictor step.
7. Corrects back to the curve using Newton's method.
8. Computes eigenvalues of the state Jacobian.
9. Detects a sign change relative to the target eigenvalue.
10. Refines the crossing by bisection.
11. Prints the detected point and the corresponding eigenvalues.

## Features

- Pure Python implementation
- No external dependencies
- Quadratic polynomial map evaluation
- Analytical Jacobian construction
- Newton correction with one fixed coordinate
- Tangent computation from the nullspace of the fixed-point Jacobian
- Eigenvalue computation for:
  - one-dimensional systems,
  - two-dimensional systems,
  - higher-dimensional systems using QR iteration
- Bisection refinement of the detected crossing

## Project structure

```text
bifurcation-continuation-solver/
├── README.md
└── bif.py
```

## Input format

The script reads all input from `stdin`.

The expected data layout is:

```text
p0_1 p0_2 ... p0_(n+1)
coefficients for f_1
coefficients for f_2
...
coefficients for f_n
direction
target_lambda
```

The script infers `n` automatically.

For a system of dimension `n`, the polynomial is evaluated in `m = n + 1` variables.  
The quadratic basis contains:

```text
1
p_1, ..., p_m
p_1 p_1, p_1 p_2, ..., p_m p_m
```

The number of coefficients for each component is:

```text
k = (n + 3)(n + 2) / 2
```

Each map component `f_i` must therefore provide exactly `k` coefficients.

## Output format

The script prints two lines:

1. The detected point on the fixed-point curve:

```text
x_1 x_2 ... x_n parameter
```

2. The eigenvalues of the state Jacobian at that point:

```text
lambda_1 lambda_2 ... lambda_n
```

All values are printed in high-precision scientific notation.

## Running the script

```bash
python bif.py < input.txt
```

Example:

```bash
python bif.py < examples/problem_01.txt
```

## Numerical method

### Predictor

A tangent vector is computed from the Jacobian of the fixed-point equation.  
The algorithm chooses a coordinate chart by fixing the coordinate with the largest tangent component.

### Corrector

Newton's method corrects the predicted point back to the fixed-point curve while holding one coordinate fixed.

### Detection

At each continuation step, the script compares the closest eigenvalue to the target value.  
When the sign changes, the crossing is refined by repeated bisection.

## Important implementation details

- The continuation step is fixed at `delta = 0.005`.
- Newton's method runs for at most 60 iterations.
- The main continuation loop allows up to 200,000 predictor-corrector steps.
- Bisection uses 100 refinement iterations.
- The target eigenvalue is parsed as an integer in the current version.

## Limitations

- The script is designed for real-valued computations.
- Complex eigenvalues are not handled explicitly.
- The QR eigenvalue routine is a simple educational implementation.
- There is limited input validation.
- The fixed step size may be inefficient for stiff or sharply curved branches.
- The script stops after detecting the first crossing.

## Future improvements

- Add adaptive step-size control.
- Add support for complex eigenvalues.
- Add command-line arguments.
- Add structured JSON or YAML input.
- Add example input files.
- Add plotting of the continuation curve.
- Replace the custom linear algebra with NumPy for robustness.
- Add tests for polynomial basis generation, Jacobians, and Newton correction.
