function [KNg,KNgi] = K_Ng(N,g)

% FUNCTION [KNg,KNgi] = K_Ng(N,g)
%
% This function builds the set of N-tuples of a Lambda-homogeneous
% polynomial of degree g.
%
% EXAMPLE:
%
%  N = [ 2 2 3 ];
%  g = [ 2 3 2 ];
%  
%  [KNg,KNgi] = K_Ng(N,g)
%
% -------------------------------------------------------------------------
%
%             JAN DE CAIGNY 25/02/2009 (KULEUVEN/UNICAMP)
%

% # of unit-simplices in the multi-simplex \Lambda
L = length(g);

% build partial N(i)-tuples of degree g(i)
for i = 1:L,
  KNgi{i} = getKd(N(i),g(i));
end

KNg = KNgi{1};
% start building all combinations of the N(i)-tuples to obtain the N-tuples
for i = 2:L,
  for k = 1:size(KNg,1),
    for j = 1:size(KNgi{i},1),
      KNgij((k-1)*size(KNgi{i},1)+j,:) = [ KNg(k,:) KNgi{i}(j,:) ];
    end
  end
  KNg = KNgij;
  clear KNgij
end

ind = 0;
for i = 1:L,
  KNgi{i} = KNg(:,ind+1:ind+N(i));
  ind     = ind+N(i);
end

% simple check
if size(KNg,1) ~= J_Ng(N,g),
  error(['Number of created N-tuples is incorrect!'])
end