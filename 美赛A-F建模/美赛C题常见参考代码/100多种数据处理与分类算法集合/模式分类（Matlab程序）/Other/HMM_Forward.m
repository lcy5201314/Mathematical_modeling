function [Pout, Alpha] = HMM_Forward(a, b, initial_state, V)

% Find the probability of a finite state in a Markov chain using the HMM forward algorithm
%
% Inputs:
%	a					- Transition probability matrix
%	b					- Output generator matrix
%	initial_state	    - Initial state or initial alpha
%	V					- Observed output sequence
%
% Output:
%	Pout				- A probability matrix
%   Alpha               - The probability matrix through the stages

alpha 		= zeros(1,size(a,1));

if (prod(size(initial_state)) == 1),
    old_alpha	= zeros(size(alpha));
    old_alpha(initial_state) = 1;
else
    old_alpha = initial_state;
end
Alpha(1,:)  = old_alpha;

%t <- t+1
for t = 2:length(V),
   for j = 1:length(alpha),
      for i = 1:size(a,1),
         %alpha_i(t-1)*a_ij*b_ij_v(t)
      	P(i,j) = old_alpha(i)*a(i,j)*b(j,V(t));   
      end
   end
   %alpha_j(t) = sum(P)
   alpha		= sum(P);
   old_alpha 	= alpha;
   Alpha(t,:)   = alpha;
end

%P(Vt) <- alpha_0(T)
Pout = alpha;
Alpha= Alpha';