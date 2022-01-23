function out = HMM_generate(a,b,initial,seq_len)

% Generate a Markov sequence
%
% Inputs:
%	a			- Transition probability matrix
%	b			- Output generator matrix
%	initial	    - Initial state
%	seq_len	    - Required sequence length
%
% Output:
%	out		    - A Markov sequence

%Check if the data is valid
if ((sum(abs(sum(a') - 1) > 1e-10) ~= 0) | (sum(abs(sum(b') - 1) > 1e-10) ~= 0)),
   error('Input matrices are incorrect')
end

%Make the sequence
out		  = zeros(1,seq_len);
curr_state = initial;

for i = 1:seq_len,
   out(i)	  = draw(b(curr_state,:));
   next_state = draw(a(curr_state,:));
end


function index = draw(probabilities)

%Select a class based on a probability vector
N 	= 1000;
P	= cumsum(probabilities);
P	= round([0 P]*N);
I	= zeros(1,N);

for i = 1:length(probabilities),
   I(P(i)+1:P(i+1)) = i;
end

%Mix the vector
I	= I(randperm(N));

index = I(1);
