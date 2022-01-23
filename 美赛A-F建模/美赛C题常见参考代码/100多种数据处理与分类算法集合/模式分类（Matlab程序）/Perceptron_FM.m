function [test_targets, a] = Perceptron_FM(train_patterns, train_targets, test_patterns, params)

% Classify using the Perceptron algorithm but at each iteration updating the worst-classified sample
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	params       	- [Maximum number of iterations, slack]
%
% Outputs
%	test_targets	- Predicted targets
%   a               - Perceptron weights
%
% NOTE: Works for only two classes.

[max_iter, slack] = process_params(params);
rate	            = 0.1;

[c, r]          = size(train_patterns);
xi  			= ones(1,r)/r*slack;

if (length(unique(train_targets)) == 2)
    train_targets = train_targets > 0;
end
train_patterns = [train_patterns ; ones(1,r)];
train_zero     = find(train_targets == 0);

%Preprocessing
y = train_patterns;
y(:,train_zero)= -y(:,train_zero);

%Initial weights
a              = sum(y')';
n			   = length(train_targets);
iter		   = 0;

while ((sum(a'*train_patterns.*(2*train_targets-1)<0)>0) & (iter < max_iter))
   iter = iter + 1;
   %Find worst-classified sample
   A            = a'*train_patterns.*(2*train_targets-1)+xi;
   [m, indice]  = min(A);
   if (a' *  y(:,indice) <= 0)
      a = a + y(:,indice);
   end
   
   %Calculate the new slack vector
   xi(indice)   = xi(indice) + rate;
   xi	        = xi / sum(xi) * slack;
   
end

if (iter == max_iter),
   disp(['Maximum iteration (' num2str(max_iter) ') reached']);
end

%Classify test patterns
test_targets = a'*[test_patterns; ones(1, size(test_patterns,2))] > 0;
a = a';