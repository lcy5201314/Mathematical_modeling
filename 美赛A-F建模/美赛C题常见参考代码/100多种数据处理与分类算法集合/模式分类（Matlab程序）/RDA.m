function test_targets = RDA (train_patterns, train_targets, test_patterns, lamda)

% Classify using the Regularized descriminant analysis (Friedman shrinkage algorithm)
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	lamda		    - Parameter for the algorithm
%
% Outputs
%	test_targets	- Predicted targets
%
% NOTE: Works for only two classes.

Uc  = unique(train_targets);
Dim = size(train_patterns, 1);

for i = 1:length(Uc),
    in = find(train_targets == Uc(i));
    
    %Estimate MLE mean and covariance
    mu      = mean(train_patterns(:,in)');
    sigma_c = cov(train_patterns(:,in)',1);
    nc      = length(in);

    %Shrink
    S       = nc * sigma_c;
    n		= nc;
    sigma   = zeros(Dim);
    nk		= n;
    sk	    = S;
   
    for j = 1:n,
       sk		 = (1 - lamda)*sk + lamda*S;
       nk		 = (1 - lamda)*nk + lamda*n;
       sigma_c   = sk / nk;
       sigma_c   = (1 - lamda) * sigma_c + lamda/2*trace(sigma_c)*eye(Dim);
       sk		 = sigma_c * nk;
   end
   
    param_struct(i).mu      = mu;
    param_struct(i).sigma   = sigma_c;
    param_struct(i).w       = 1;
    param_struct(i).p       = 1/length(Uc);
    param_struct(i).type    = 'Gaussian';
end
   
test_targets = classify_paramteric(param_struct, test_patterns);

