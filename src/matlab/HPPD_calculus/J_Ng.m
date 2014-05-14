function JNg = J_Ng(N,g)

% FUNCTION JNg = J_Ng(N,g)
%
% This function calculates the number of N-tuples of a Lambda-homogeneous
% polynomial of degree g.
%
% EXAMPLE:
%
%  N = [ 2 2 3 ];
%  g = [ 2 3 2 ];
%  
%  JNg = J_Ng(N,g)
%
% -------------------------------------------------------------------------
%
%             JAN DE CAIGNY 25/02/2009 (KULEUVEN/UNICAMP)
%

JNg = 1;
for i = 1:length(g),
  JNgi = getJd(N(i),g(i));
  JNg  = JNg * JNgi;
end