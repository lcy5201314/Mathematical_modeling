function test_targets = parzen(train_patterns, train_targets, test_patterns, hn)

% Classify using the Parzen windows algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	hn      	    - Normalizing factor for h
%
% Outputs
%	test_targets	- Predicted targets

N       = size(test_patterns, 2);
Uc      = unique(train_targets);
V		= zeros(length(Uc), N);
x_i     = train_patterns;

for j = 1:length(Uc),
    indices = find(train_targets == Uc(j));
    P(j)    = length(indices)/size(x_i,2);
    n		= length(indices);
    
    for i = 1:n,
        temp      = sum((test_patterns - x_i(:,indices(i))*ones(1,N)).^2);
        V(j,:)    = V(j,:) + phi(temp./hn);
    end
    
    V(j,:) = V(j,:) / sum(V(j,:)) * P(j);
end

test_targets = zeros(1, N);
for i = 1:N,
    [m, best]       = max(V(:,i));
    test_targets(i) = Uc(best);
end
%END Parzen

function p = phi(val)

%The window function for the Parzen window
p = (abs(val) <= 0.5);