function test_targets = ML_II(train_patterns, train_targets, test_patterns, Ngaussians)

% Classify using the ML-II algorithm. This function accepts as inputs the maximum number
% of Gaussians per class and returns a decision surface based on the most likely number of
% Gaussians in each class
%
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	Ngaussians      - Maximum number of Gaussians per class
%
% Outputs
%	test_targets	- Predicted targets
%

Uc          = unique(train_targets);
[Ndim, M]	= size(train_patterns);
P_D_given_h	= zeros(length(Uc), max(Ngaussians));

%Use holdout for diffrentiating between training data for finding the gaussian parameters
%and the likelihood
holdout			= 0.5;
i				= randperm(M);
train_patterns  = train_patterns(:,i);
train_targets	= train_targets(i);
EMindices		= 1:floor(M*holdout);
MLindices		= floor(M*holdout)+1:M;
i0				= MLindices(find(train_targets(MLindices) == 0));
i1				= MLindices(find(train_targets(MLindices) == 1));
Ni0				= length(i0);
Ni1				= length(i1);

for i = 1:max(Ngaussians),
    
    %Calculate the error for each possible model. Use only EMindices data
    G			= min([ones(1,length(Ngaussians))*i; Ngaussians]);
    [D, param]	= EM(train_patterns(:,EMindices), train_targets(:,EMindices), train_patterns(:,1:2), G);
    
    %Calculate likelihood of the data given these Gaussians
    %Use only the MLindices data
    for j = 1:length(Uc)
        if (P_D_given_h(j, G(j)) == 0),	%Do it only if it wasn't already computed
            in = MLindices(find(train_targets(MLindices) == Uc(j)));
            P_D_given_h(j, G(j)) = computeML(train_patterns(:,in), param(j).mu, param(j).sigma, param(j).w);      
        end
    end
end

[best, bestNgaussians]  = min(P_D_given_h');
test_targets            = EM(train_patterns, train_targets, test_patterns, bestNgaussians);



function P = computeML(patterns, mu, sigma, w)

M	= size(patterns,2);
Ng	= size(mu,1);
p	= zeros(Ng, M);

warning off

if Ng == 1,
    for j = 1:M,
        x		= patterns(:,j);
        p(j)	= w(1)/(2*pi*sqrt(det(sigma)))*exp(-0.5*(x-mu')'*inv(sigma)*(x-mu'));
    end
    P			= prod(p);
else
    for j = 1:M,
        x		= patterns(:,j);
        for k = 1:length(w),
            p(k, j) = w(k)/(2*pi*sqrt(det(squeeze(sigma(k,:,:)))))*...
                exp(-0.5*(x-mu(k,:)')'*inv(squeeze(sigma(k,:,:)))*(x-mu(k,:)'));
        end
    end
    P			= prod(sum(p));   
end

