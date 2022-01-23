function test_targets = Minimum_Cost(train_patterns, train_targets, test_patterns, lambda)

% Classify using the minimum error criterion via histogram estimation of the densities
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	lambda  - Cost matrix
%
% Outputs
%	test_targets	- Predicted targets
%
% Note: Works for only two classes

Nbins             = max(3,floor(size(train_patterns,2).^(1/3)));
Uc                = unique(train_targets);
[H, Bins, region] = high_histogram(train_patterns, Nbins);

for i = 1:length(Uc),
    in   = find(train_targets == Uc(i));
    P(i) = length(in) / length(train_targets);
    p(i,:,:) = high_histogram(train_patterns(:,in),Nbins,region);
end

decision    = (lambda(2,1) - lambda(1,1))*P(1)*squeeze(p(1,:,:)) < ...
              (lambda(1,2) - lambda(2,2))*P(2)*squeeze(p(2,:,:));

%Classify test samples
[H, Bins, region] = high_histogram(test_patterns, Nbins, region);
test_targets      = zeros(1, size(Bins,2));
for i = 1:size(Bins,2),
    test_targets(i) = decision(Bins(1,i),Bins(2,i));
end
