function Path = HMM_Decoding(a, b, initial_state, V)

% Find the most likely states from the outputs of a Markov chain using the HMM decoding algorithm
%
% Inputs:
%	a					- Transition probability matrix
%	b					- Output generator matrix
%	initial_state	    - Initial state
%	V					- Observed output sequence
%
% Output:
%	Path				- The most likely sequences

T		= length(V);
c		= size(a,1);
Path	= initial_state;

alpha   = zeros(1,c);
alpha(initial_state) = 1;

%t <- t + 1
for t = 2:T,
    old_alpha   = alpha;
    alpha	    = zeros(1,c);
    
    %alpha_k(t) <- b_jkv(t)*sum(alpha_i(t-1)*a_ij    
    for k = 1:c,
        for i = 1:c,
            P(i) = old_alpha(i)*a(i,k);   
        end
        alpha(k) = b(k,V(t)) * sum(P);
    end
    
    %j' <- argmax(alpha)
    [m, j] = max(alpha);
    
    %Append to path
    Path = [Path j];
end
