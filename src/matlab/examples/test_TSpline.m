
clear all
close all
clear global
clc

addpath ../Bspline


alpha = 1;      
knots = [0, 0, pi/2, pi, pi];
degree = 2;
x = linspace(0,pi,1e3)';

i1 = (x <= pi/2);
i2 = (x >= pi/2);
B1 = zeros(size(x));
B1(i1) = sin(pi/2-x(i1));
B2 = zeros(size(x));
B2(i1) = sin(x(i1));
B2(i2) = sin(pi-x(i2));
B3 = zeros(size(x));
B3(i2) = sin(x(i2)-pi/2);
B1 = [B1,B2,B3];

b = TSplineBasis(2*knots, degree);
B2 = b.f(2*x);

figure
subplot(121)
plot(x,B1)
subplot(122)
plot(x,B2)

