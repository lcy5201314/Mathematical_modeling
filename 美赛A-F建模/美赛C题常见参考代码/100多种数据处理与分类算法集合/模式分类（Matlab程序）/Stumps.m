function [test_targets, w] = Stumps(train_patterns, train_targets, test_patterns, params)

% Classify using simple stumps algorithm
% Inputs:
% 	train_patterns  - Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test patterns
%	weights	        - Unused (Except if weighted stumps is needed)
%
% Outputs
%	test_targets    - Predicted labels 
%	w			    - Decision surface parameters
%
% NOTE: Works for only two classes!

train_one  = find(train_targets == max(train_targets));
train_zero = find(train_targets == min(train_targets));

if (length(params)-1 == length(train_targets)), 
    p = params(1:end-1);
else
    p = ones(size(train_targets));   
end

dim        = size(train_patterns,1);
w          = zeros(1,dim);
err        = zeros(1,dim);
direction  = zeros(1,dim);

for i = 1:dim,
    %For each dimension, find the point where a stump gives the minimal error
    
    %First, sort the working dimension 
    [data(i,:), indices] = sort(train_patterns(i,:));
    temp_targets         = train_targets(indices);
    temp_p		 	     = p(indices);
    
    decision             = cumsum(temp_p .* temp_targets)/length(train_one) - cumsum(temp_p .* (~temp_targets))/length(train_zero);
    [err(i),W]           = max(abs(decision));
    w(i)                 = data(i,W);
    direction(i)         = sign(decision(W));
end

[m, min_dim] = max(err);
indices      = find(~ismember(1:dim,min_dim));
w(indices)   = 0;

if (direction(min_dim) > 0)
    indices = find(test_patterns(min_dim,:) < w(min_dim));
else
    indices = find(test_patterns(min_dim,:) > w(min_dim));
end

test_targets          = zeros(1, size(test_patterns,2));
test_targets(indices) = 1;

