function [test_targets, g0, g1] = NDDF(train_patterns, train_targets, test_patterns, cost)

% Classify using the normal density discriminant function
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	cost			- Cost for class 0 (Optional, Unused yet)
%
% Outputs
%	test_targets	- Predicted targets
%	g    	        - The discriminant function for the test examples

[d, L] = size(train_patterns);

%Estimate mean and covariance for each class
mu		= zeros(d,length(unique(train_targets)));
sigma	= zeros(d,d,length(unique(train_targets)));
p		= zeros(length(unique(train_targets)));

classes = unique(train_targets);
for i = 1:length(classes),
   indices			= find(train_targets == classes(i));
   mu(:,i)			= mean(train_patterns(:,indices)')';
   sigma(:,:,i)	    = cov(train_patterns(:,indices)',1)';
   p(i)				= length(indices)/length(train_targets);
end

test_targets = zeros(1, size(test_patterns,2));
for i = 1:size(test_patterns,2),
    sample = test_patterns(:,i);
    for j = 1:length(classes),
        g(i,j) = -0.5*(sample - mu(:,j))'*inv(squeeze(sigma(:,:,j)))*(sample - mu(:,j)) - ...
   		          d/2*log(2*pi)-0.5*log(det(squeeze(sigma(:,:,j))))+log(p(j));
    end
    [m, best]       = max(g(i,:));
    test_targets(i) = classes(best);
end