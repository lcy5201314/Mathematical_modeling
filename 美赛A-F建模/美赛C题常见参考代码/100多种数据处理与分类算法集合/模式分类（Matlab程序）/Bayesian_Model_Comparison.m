function test_targets = Bayesian_Model_Comparison(train_patterns, train_targets, test_patterns, Ngaussians)

% Classify using the Bayesian model comparison algorithm. This function accepts as inputs
% the maximum number of Gaussians per class and returns a decision surface based on the 
% most likely number of Gaussians in each class
%
% Inputs:
% 	training_patterns   - Train patterns
%	training_targets	- Train targets
%   test_patterns       - Test  patterns
%	Ngaussians          - Number of redraws
%
% Outputs
%	test_targets        - Predicted targets
%
% Strongly built for only two classes!

[Ndim, M]	= size(train_patterns);
P_D_given_h	= zeros(2, max(Ngaussians)+2);

%Use holdout for diffrentiating between training data for finding the gaussian parameters
%and the likelihood
holdout			    = 0.5;
i					= randperm(M);
train_patterns      = train_patterns(:,i);
train_targets	    = train_targets(i);
EMindices		    = 1:floor(M*holdout);
MLindices		    = floor(M*holdout)+1:M;
i0					= MLindices(find(train_targets(MLindices) == 0));
i1					= MLindices(find(train_targets(MLindices) == 1));
Ni0				    = length(i0);
Ni1				    = length(i1);

for i = 1:max(Ngaussians),
   
   %Calculate the error for each possible model. Use only EMindices data
   G			= min([ones(1,length(Ngaussians))*i; Ngaussians]);
   [D, param]	= EM(train_patterns(:,EMindices), train_targets(:,EMindices), train_patterns(:,EMindices), G);
   
   %Calculate likelihood of the data given these Gaussians
   %Use only the MLindices data
   if (P_D_given_h(1, G(1)) == 0),	%Do it only if it wasn't already computed
		P_D_given_h(1, G(1)) = computeML(train_patterns(:,i0), param(1).mu, param(1).sigma, param(1).w);      
   end
   if (P_D_given_h(2, G(2)) == 0),	%Do it only if it wasn't already computed
		P_D_given_h(2, G(2)) = computeML(train_patterns(:,i1), param(2).mu, param(2).sigma, param(2).w);      
   end
   
end
P_D_given_h(find(isnan(P_D_given_h)))	=	0;
P_D_given_h = P_D_given_h./(eps+sum(P_D_given_h')'*ones(1,size(P_D_given_h,2))); %Normalize

%Compute the Hessian for each class
H 						= diff([zeros(2,2), P_D_given_h], 2, 2);

P_D_given_h			= P_D_given_h .* (abs(H).^(-0.5)) .* (2*pi);

likelihood = P_D_given_h(1,[1:Ngaussians(1)])' * P_D_given_h(2,[1:Ngaussians(2)]);

%Choose the ML model as the one with the lowest error
[i1, i2] = find(likelihood == max(max(likelihood)));
if isempty(i1),
    error('Could not find a likely pair.')
end

i1       = i1(1); i2 = i2(1); %To give preference for simpler models...

test_targets = EM(train_patterns, train_targets, test_targets, [i1, i2]);

disp(['FINAL SELECTION: Using ' num2str(i1) ' Gaussians for class 1 and ' num2str(i2) ' Gaussians for class 2'])


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

