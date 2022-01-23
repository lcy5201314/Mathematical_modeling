function [test_targets, E] = ada_boost(train_patterns, train_targets, test_patterns, params)

% Classify using the AdaBoost algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	Params	- [NumberOfIterations, Weak Learner Type, Learner's parameters]
%
% Outputs
%	test_targets	- Predicted targets
%   E               - Errors through the iterations
%
% NOTE: Suitable for only two classes
%

[k_max, weak_learner, alg_param] = process_params(params);

[Ni,M]			= size(train_patterns);
W			 	= ones(1,M)/M;
IterDisp		= 10;

full_patterns   = [train_patterns, test_patterns];
test_targets    = zeros(1, size(test_patterns,2));

%Do the AdaBoosting
for k = 1:k_max,
   %Train weak learner Ck using the data sampled according to W:
   %...so sample the data according to W
   randnum = rand(1,M);
   cW	   = cumsum(W);
   indices = zeros(1,M);
   for i = 1:M,
      %Find which bin the random number falls into
      loc = max(find(randnum(i) > cW))+1;
      if isempty(loc)
         indices(i) = 1;
      else
         indices(i) = loc;
      end
   end
   
   %...and now train the classifier
   Ck 	= feval(weak_learner, train_patterns(:, indices), train_targets(indices), full_patterns, alg_param);

   %Ek <- Training error of Ck 
   E(k) = sum(W.*(Ck(1:M) ~= train_targets));
   
   if (E(k) == 0),
      break
   end
   
   %alpha_k <- 1/2*ln(1-Ek)/Ek)
   alpha_k = 0.5*log((1-E(k))/E(k));
   
   %W_k+1 = W_k/Z*exp(+/-alpha)
   W  = W.*exp(alpha_k*(xor(Ck(1:M),train_targets)*2-1));
   W  = W./sum(W);
   
   %Update the test targets
   test_targets  = test_targets + alpha_k*(2*Ck(M+1:end)-1);
   
   if (k/IterDisp == floor(k/IterDisp)),
      disp(['Completed ' num2str(k) ' boosting iterations'])
   end
   
end

test_targets = test_targets > 0;