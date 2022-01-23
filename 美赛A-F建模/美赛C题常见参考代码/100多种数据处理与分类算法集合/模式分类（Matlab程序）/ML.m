function test_targets = ML(train_patterns, train_targets, test_patterns, AlgorithmParameters)

% Classify using the maximum-likelyhood algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	params  		- Unused
%
% Outputs
%	test_targets	- Predicted targets

Uclasses = unique(train_targets);

for i = 1:length(Uclasses),
    indices = find(train_targets == Uclasses(i));
    
    %Estimate mean and covariance
    param_struct(i).mu      = mean(train_patterns(:,indices)');
    param_struct(i).sigma   = cov(train_patterns(:,indices)',1);
    param_struct(i).p       = length(indices)/length(train_targets);
    param_struct(i).w       = 1/length(Uclasses);
    param_struct(i).type    = 'Gaussian';
end

%Classify test patterns
test_targets = classify_paramteric(param_struct, test_patterns);