function [Pout, Beta] = HMM_Backward(a, b, final_state, V)

% Find the probability of a finite state in a Markov chain using the HMM backward algorithm
%
% Inputs:
%	a					- Transition probability matrix
%	b					- Output generator matrix
%	final_state		    - Final state or final beta
%	V					- Observed output sequence
%
% Output:
%	Pout				- A probability matrix
%   Beta                - The probability matrix through the stages


beta 		= zeros(1,size(a,1));

if (prod(size(final_state)) == 1),
    old_beta	= zeros(size(beta));
    old_beta(final_state) = 1;
else
    old_beta = final_state;
end

Beta(length(V),:)   = old_beta;

%t <- t-1
for t = length(V)-1:-1:1,
   for j = 1:length(beta),
      for i = 1:size(a,1),
         %beta_i(t+1)*a_ij*b_ij_v(t)
      	P(i,j) = old_beta(i)*a(i,j)*b(j,V(t+1));   
      end
   end
   %beta_j(t) = sum(P)
   beta			= sum(P);
   old_beta 	= beta;
   Beta(t,:)= beta;
end

%P(Vt) <- beta_i(0)
Pout = beta;
Beta = Beta';
