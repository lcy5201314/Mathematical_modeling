function P = HMM_evaluation(a, b, V)

% Find the probability of a finite state in a Markov chain 
%
% Inputs:
%	a					- Transition probability matrix
%	b					- Output generator matrix
%	V					- Observed output sequence
%
% Output:
%	P					- Probability of the sequence

lenV = length(V);

%First, build a matrix with all the possible sequences of length(len(V))
Sequence = build_sequence(lenV, size(a,1));

p			= zeros(1,size(Sequence,1));
for i = 1:size(Sequence,1),
   for j = 1:size(Sequence,2),
      if (j == 1),
         p(i) = b(Sequence(i,1),V(1));
      else
         p(i) = p(i)*a(Sequence(i,j),Sequence(i,j-1))*b(Sequence(i,j-1),V(j));
      end
   end
end

P = sum(p);



function S = build_sequence(depth, Nstates)
%Build all possible sequences recursively. 
%depth is the remaining depth and Nstates are the number of possible sattes in the chain
S = [];

if (depth > 1),
   S1 = build_sequence(depth-1, Nstates);
   for i = 1:Nstates,
      S = [S; i*ones(size(S1,1),1) S1];
   end
else
   S = [1:Nstates]';
end

         
         