# Splines
A framework for manipulating and optimizing (B-)splines

## Installation
Simply place the `src` folder in your MATLAB path

### Dependencies
For optimizing splines [YALMIP](http://users.isy.liu.se/johanl/yalmip/) is
required. For the more adventurous, there is basic support for
[CasADi](http://www.casadi.org), which is much more suited for nonconvex optimization
problems.

## Basic usage
The code below shows a basic example how a parametric optimization problem can
be solved approximately with Splines.

```#!matlab
d = 3;   % degree
L = 4;   % Range [0, L]
n = 21;  % Number of knots
Bl = BSplineBasis([0 L], d, n);
l = parameter();
x = BSpline.sdpvar(Bl, [2, 1]);  % x is 2x1 spline function with basis Bl
obj = x(1) + 2 * x(2);
con = [x(1) >= 0, x(2) >= 0, x(2) <= 2, x(1) + l * x(2) <= 2];
sol = optimize(con, -obj.integral);
obj = value(obj);
```

Those familiar with YALMIP should notice the strong similarity with its syntax.

### Basis
A B-spline basis is defined by a non-decreasing knot sequence `knots` and a degree
`deg`. To construct such a basis, use

    >> B = BSplineBasis(knots, deg);

A uniform B-spline basis of degree `deg` with `n` breakpoint on an interval
`I` can be created with

    >> B = BSplineBasis(I, deg, n);

### Coefficients
The coefficients for a spline are created by the `Coefficients` class

    >> coeffs = Coefficients(data, siz, shape);

where `data` contains the numerical data of the coefficients, `siz` indicates
the size of the coefficients tensor and `shape` is the size of an individual
coefficient. For example, the following code creates a 3x3x3 coefficient tensor
of 2x2 coefficients.

    >> data = reshape(1:108, [6, 6, 3]);
    >> coeffs = Coefficients(data, [3, 3, 3], [2, 2]);
    >> coeffs(1).data

    ans =

         1     7
         2     8

### Constructing a B-spline
A (tensor product) spline is defined by a basis for each dimension and a
coefficient tensor. Furthermore a B-spline can be scalar valued or matrix
valued. A spline is created using

    >> S = BSpline(bases, coeffs);

where bases is a 1 x n cell-array of `BSplineBasis` instances and coeffs is
either an instance of `Coefficients` or a cell array of coefficients of m_1 x
m_2 x ... x m_n where m_i is the length of the corresponding `BSplineBasis`.
As this definition is rather cumbersome to implement, two shortcuts are
available.

1) For scalar valued splines the coefficients can be supplied as a regular m_1
   x m_2 x ... x m_n matrix instead of a cell array.

2) Furthermore, for univariate splines the bases is allowed to be a
   `BSplineBasis` object in stead of a cell array.

Examples:

A univariate spline:

    >> B1 = BSplineBasis([0, 0, 0, 0.5, 1, 1, 1], 2);
    >> c1 = randn(length(B1), 1);
    >> S1 = BSpline(B1, c1);

A multivariate scalar valued spline:

    >> B2 = BSplineBasis([0, 0, 0, 0, 0.33, 0.66, 1, 1, 1, 1], 3);
    >> c2 = randn(length(B1), length(B2));
    >> S2 = BSpline({B1, B2}, c2);

### Operations on splines
For the various methods that can be called on splines, refer to the `examples`
folder in the repository or have a look at the `Function.m` and `BSpline.m`
files in the `src` folder.

### Optimizing over splines
To create a spline optimization variable use

    >> S = BSpline.sdpvar(bases, shape[, props]);

where `bases` contains a cell array of the B-spline basis that construct the
tensor product basis, `shape` is the dimension of the individual coefficients
and `props` are additional YALMIP properties assigned to the coefficients,
e.g. `'symmetric'` for symmetrical matrix valued coefficients.

Constraints on splines are guaranteed by imposing simple constraints on the
coefficients. Constraints are _not_ sampled. Constraints on splines are
imposed using the relational operators `<=` and `>=`.

    >> con = S >= 0;

which returns a list of YALMIP constraints.
