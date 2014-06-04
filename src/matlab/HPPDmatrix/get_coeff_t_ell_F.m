function c = get_coeff_t_ell_F(g,F)

% FUNCTION c = coef_t_ell_F_unit(g,F)
%
% Given the linear relationship \alpha = F \gamma, where F is a matrix of
% dimension N x M, between the variable \alpha that takes values in the
% unit-simplex of dimension N and the variable \gamma that takes values in
% the unit-simplex of dimension M, this file calculates the coefficients
% that are necessary to calculate the linear combination for the change of
% variables between the coefficients of the homogeneous polynomials
% A(\alpha) = Ah(\gamma) of degree g.
%
% USE: c = coef_t_ell_F_unit(g,F)
%
%      g : degree of the homogeneous polynomial
%      F : matrix that expresses the linear relation between \alpha and
%          \gamma
%            size : N x M
%
%      c : coefficients for the linear combination
%            size : J_M(g) x J_N(g)
%
% EXAMPLE: 
%
%  N = 3;
%  g = 2;
%  b = 0.2;
%  % get linear relation between \alpha and \gamma
%  [F,H,FpH,M] = buildGamma(N,b);
%
%  c = coef_t_ell_F_unit(g,F{1})
%
% -------------------------------------------------------------------------
%
%             JAN DE CAIGNY 25/06/2009 (KULEUVEN/UNICAMP)
%

%%%%%%%%%%%%%%%%%
% PREPROCESSING %
%%%%%%%%%%%%%%%%%

  % size of the unit-simplices
  [N,M] = size(F);
  
  % create the vector f that consists of all rows of F
  f = vec(F')';
  
  % create dimension of the multi-simplex (M,M,...,M) \in \N^N
  NM = M*ones(1,N);

  % # of monomials in the alpha representation
  JNg = get_JNg(N,g);
  
  % # of monomials in the gamma representation
  JMg = get_JNg(M,g);

  % corresponding monomial coefficients
  KNg = get_KNg(N,g);
  KMg = get_KNg(M,g);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD NEW REPRESENTATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % LOOP 1: over t, the power of gamma
  for t = 1:JMg, tt = KMg(t,:);    
    % LOOP 2: over l, the power of alpha
    for l = 1:JNg, ll = KNg(l,:);
      
      % define the current coefficient as zero
      c(t,l) = 0;
      
      % get multi-simplex monomial powers for gamma
      JNMll          = J_Ng(NM,ll);
      [KNMll,KNMllj] = K_Ng(NM,ll);
      % LOOP 3: inner loop over all possible values of k
      for k = 1:JNMll, kk = KNMll(k,:);
        % get sum of k_i
        sum_kj = zeros(1,M);
        for j = 1:N
          sum_kj = sum_kj + KNMllj{j}(k,:);
        end
        % compare tt with sum_k1
        if sum(tt == sum_kj) == length(tt),
          c(t,l) = c(t,l) + Pi(ll)/Pi(kk) * prod(f.^kk);
        end
      end
    end
  end