# Explicitely parameterized approximate solutions of parametric cone programs

This software provides MATLAB and Python tools for solving parametric cone
programs using B-splines.

For example the following code solves a simple parametric LP

```
#!matlab
deg = 3;
L = 4;
n = 21;
Bl = BSplineBasis([0 * ones(1, deg), linspace(0, L, n), L * ones(1, d)], deg);
l = BSpline(BSplineBasis([0, 0, L, L], 1), [0, L]');

c1 = sdpvar(length(Bl), 1);
c2 = sdpvar(length(Bl), 1);
x1 = BSpline(Bl, c1);
x2 = BSpline(Bl, c2);
obj = x1 + 2 * x2;
options = sdpsettings('verbose',1);
sol = solvesdp([x1 >= 0, x2 >= 0, x2 <= 2, x1 + l * x2 <= 2], -obj.integral, options);
```

A working document on the theory is found in the doc folder of this repository.

## Matlab docs

The documentation shows the basic syntax for defining BSplines.

### Installation

Add the folder Function and its subfolders to the MATLAB path and you are
ready to go.

### Constructing a basis

A polynomial basis is defined simply by its degree `deg` and is construct with

    Q = PolynomialBasis(deg);

A B-spline basis is defined by a non-decreasing knot sequence `knots` and a degree
`deg`. To construct such a basis, use

    B = BSplineBasis(knots, deg);

### Constructing a B-spline

A (multivariate) B-spline is defined by a basis for each variable and a
coefficient tensor. Furthermore a B-spline can be scalar valued or matrix
valued. In MATLAB, a B-Spline is created by the syntax

   S = BSpline(bases, coefficients)

where bases is a 1 x n cell-array of `BSplineBasis` objects and coefficients
is a cell array of coefficients of m_1 x m_2 x ... x m_n where m_i is the
length of the corresponding `BSplineBasis`. As this definition is rather
cumbersome to implement, two shortcuts are available.

1) For scalar valued splines the coefficients can be supplied as a regular m_1
   x m_2 x ... x m_n matrix in stead of a cell array.

2) Furthermore, for univariate splines the bases is allowed to be a
   `BSplineBasis` object in stead of a cell array.

Examples:

A univariate spline:

    B1 = BSplineBasis([0, 0, 0, 0.5, 1, 1, 1], 2);
    c1 = randn(length(B1), 1);
    S1 = BSpline(B1, c1);

A multivariate scalar valued spline:
    
    B2 = BSplineBasis([0, 0, 0, 0, 0.33, 0.66, 1, 1, 1, 1], 3);
    c2 = randn(length(B1), length(B2));
    S2 = BSpline({B1, B2}, c2);

A multivariate matrix valued 3x3 spline:

    c3 = cell(length(B1), length(B2));
    for i=1:length(B1)
        for j=1:length(B2)
            c3{i, j} = randn(3, 3);
        end
    end
    S3 = BSpline({B1, B2}, c3);

### Manipulating B-Splines

Summation, multiplication, concatenation and (conjugated) transposition are
overloaded in the BSpline class. Therefore, the regular MATLAB syntax `S1 +
S2`, `S1 * S2`, `[S1, S2]` and `S1'` can be used. Furthermore, the instance
method `integrate` returns the definite integral over the spline domain and
`derivative` returns the derivative of the B-spline as a `BSpline` object.

### Constructing a polynomial

Polynomials are defined using a similar syntax to B-splines:

    P = Polynomial(bases, coefficients);