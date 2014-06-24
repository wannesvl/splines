
clear all
close all
clear global
clc

addpath ../Function/
addpath ../Function/Coefficients/
addpath ../Function/Basis


V = [-1, -1, 0, 1,  1,  0 ;         % [x;y] of vertices
      0,  1, 1, 0, -1, -1 ];
  
N = 4;      % number of triangles 
Bb = BSplineBasis([0, 0, 1, 1], 1);
Ba = BSplineBasis([0, linspace(0,1,N+1), 1], 1);

Cx = [ V(1,1), V(1,1), V(1,6), V(1,6), V(1,5) ;
       V(1,2), V(1,3), V(1,3), V(1,4), V(1,4) ];
x = BSpline({Bb,Ba}, Cx);   

Cy = [ V(2,1), V(2,1), V(2,6), V(2,6), V(2,5) ;
       V(2,2), V(2,3), V(2,3), V(2,4), V(2,4) ];
y = BSpline({Bb,Ba}, Cy);

n = 101;
a = linspace(0, 1/N, n)'*ones(1,N) + 1/N*ones(n,1)*(0:N-1);
b = linspace(0,1,n)';

figure
for i = 1:N
    xf = x.f({b,a(:,i)});
    yf = y.f({b,a(:,i)});
    plot(xf(:), yf(:), '.'), hold all
end
plot([V(1,:),V(1,1)], [V(2,:),V(2,1)], 'k')