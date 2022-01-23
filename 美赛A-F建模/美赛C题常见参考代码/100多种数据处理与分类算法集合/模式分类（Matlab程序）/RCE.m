function test_targets = RCE(train_patterns, train_targets, test_patterns, lambda_m)

% Classify using the reduced coulomb energy algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	lambda_m	    - Maximum radius 
%
% Outputs
%	test_targets	- Predicted targets

epsilon 	= 1e-4;
[Dim,Nf]    = size(train_patterns);
X           = train_patterns;

%Train the classifier
W   		= X;            %w_ij <- x_i
lambda      = zeros(1,Nf);

for i = 1:Nf,
    %x_hat <- arg min D(x, x_tag)
    dist 			 = sqrt(sum((X - X(:,i) * ones(1,Nf)).^2));
    [m, indices]      = sort(dist);
    x_hat			 = find(train_targets(indices) ~= train_targets(i));
    
    %lambda_j <- min(D(x_hat, x_tag)-epsilon, lambda_m)
    lambda(i)	     = min(dist(x_hat(1))-epsilon,lambda_m);
end

%Classify the test samples
N               = size(test_patterns, 2);
Uc              = unique(train_targets);
test_targets    = zeros(1, N);

for i = 1:N,
    %If D(x, x_hat_j)<lambda_j then D_t <- D_t U x_tag_j
    dist	= sqrt(sum((X - test_patterns(:,i) * ones(1,Nf)).^2));
    indices = find(dist < lambda);

    %The decision is a little different from DH&S, since there can be an ambiguous result.
    %Here we do not allow this.
    if isempty(indices),
        test_targets(i) = rand(1) > .5;
    else
        vote            = hist(train_targets(indices), Uc);
        [m, best]       = max(vote);
        test_targets(i) = Uc(best);        
    end    
end