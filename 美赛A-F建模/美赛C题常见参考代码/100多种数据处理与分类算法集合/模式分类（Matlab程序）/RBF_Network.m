function [test_targets, mu, Wo] = RBF_Network(train_patterns, train_targets, test_patterns, Nh)

% Classify using a radial basis function network algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	Nh              - Number of hidden units
%
% Outputs
%	test_targets	- Predicted targets
%   mu              - Hidden unit locations
%   Wo              - Output unit weights

[Ni, M] = size(train_patterns);

%First, find locations for the hidden unit centers using k-means
Npoints                     = 100;
[mu, center_targets, label] = k_means(train_patterns, train_targets, Nh, 0);

%Remove bad centers
ok  = find(isfinite(mean(mu)));
mu  = mu(:, ok);
Nh  = length(ok);

%Variance of the gaussians
dist = zeros(Nh);
for i=1:Nh,
    dist(i,:) = sqrt(sum((mu(:,i)*ones(1,Nh) - mu).^2));
end
max_dist = max(max(dist));
sigma    = max_dist/sqrt(2*Nh);

%Compute the activation for each pattern at each center
Phi = zeros(Nh, M);
for i = 1:Nh,
    Phi(i,:) = 1/(2*pi*sigma^2)^(Ni/2)*exp(-sum((train_patterns-mu(:,i)*ones(1,M)).^2)/(2*sigma^2));
end

%Now, find the hidden to output weights
Wo  = (pinv(Phi)'*(train_targets*2-1)')';

%Classify test patterns
N   = size(test_patterns, 2);
Phi = zeros(Nh, N);
for i = 1:Nh,
    Phi(i,:) = 1/(2*pi*sigma^2)^(Ni/2)*exp(-sum((test_patterns-mu(:,i)*ones(1,N)).^2)/(2*sigma^2));
end

test_targets = Wo * Phi;

%If there are only two classes, collapse them
if (length(unique(train_targets)) == 2)
    test_targets = test_targets > 0.5;
end