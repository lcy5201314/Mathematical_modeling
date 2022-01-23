function [test_targets, a_plus, a_minus] = Balanced_Winnow(train_patterns, train_targets, test_patterns, params)

% Classify using the balanced Winnow algorithm
% Inputs:
% 	training_patterns   - Train patterns
%	training_targets	- Train targets
%   test_patterns       - Test  patterns
%	params		        - [Num iter, Alpha, Convergence rate]
%
% Outputs
%	test_targets        - Predicted targets
%   a_plus	            - The positive weight vector
%   a_minus	            - The negative weight vector

[c, r]		    		  = size(train_patterns);
[Max_iter, alpha, eta]    = process_params(params);
y 			              = [train_patterns ; ones(1,r)];
z        			      = train_targets;

%Initial weights
a_plus	        = sum(y')';
a_minus	        = -sum(y')';
iter  	        = 0;

while (iter < Max_iter)
   iter = iter + 1;
   
   for k = 1:r,
      if (sign(a_plus'*y(:,k) - a_minus'*y(:,k)) ~= sign(z(k)-.5)),
         if (z(k) == 1),
            a_plus	= alpha.^y(:,k).*a_plus;
            a_minus	= alpha.^-y(:,k).*a_minus;
         else
            a_plus	= alpha.^-y(:,k).*a_plus;
            a_minus	= alpha.^y(:,k).*a_minus;
         end
      end
   end
   
end

a = (a_plus + a_minus)/2;

%Classify the test patterns
test_targets = a' * [test_patterns; ones(1,size(test_patterns,2))];

if (length(unique(train_targets))== 2)
    test_targets = test_targets < 0.5;
end