function b = compk(k1,k2)

% FUNCTION b = compk(k1,k2)
%
% This function compares two N-tuples of equal length and returns 1 if k1
% is element-wise bigger than or equal to k2 and 0 if not.
%
% EXAMPLE:
%
% k1 = [ 1 0 2 6 3 ];
% k2 = [ 0 0 2 5 1 ];
% k3 = [ 1 1 1 1 1 ];
%
% compk(k1,k2), compk(k1,k3)
%
% -------------------------------------------------------------------------
%
%             JAN DE CAIGNY 25/02/2009 (KULEUVEN/UNICAMP)
%

if length(k1)==length(k2)
    b=0;
    if sum(k1>=k2)==length(k1),
        b=1;
    end
else
    error('Cannot compare monomial coefficients of different length!')
end