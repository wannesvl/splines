function JNg = get_JNg(N,g)

% FUNCTION Jd = getJd(N,d)
%
%     impementation based on:
% "Stability of polytopes of matrices via affine parameter-dependent
% Lyapunov functions: Asymptotically exact LMI conditions," Ricardo C.L.F.
% oliveira, Pedro L.D. Peres, Linear Algebra and its Applications 405
% (2005), 209-228.
%
% This function calculates the number of N-tuples of a polynomial of degree
% d depending on N variables.
%
% -------------------------------------------------------------------------
%
%             Jan De Caigny 17/06/2008 (KULEUVEN/UNICAMP)
%

JNg = factorial(N+g-1)/(factorial(g)*factorial(N-1));