function [test_targets, a, updates] = Perceptron_Batch(train_patterns, train_targets, test_patterns, params)

% Classify using the batch Perceptron algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	param		    - [Num iter, Theta (Convergence criterion), Convergence rate]
%
% Outputs
%	test_targets	- Predicted targets
%   a               - Perceptron weights
%   updates         - Updates throughout the learning iterations
%
% NOTE: Works for only two classes.

[c, r]		   = size(train_patterns);
[Max_iter, theta, eta] = process_params(params);

train_patterns = [train_patterns ; ones(1,r)];
train_one      = find(train_targets == 1);
train_zero     = find(train_targets == 0);

%Preprocessing
y              = train_patterns;
y(:,train_zero)= -y(:,train_zero);

%Initial weights
a              = sum(y')';
iter  	       = 0;

update		   = 10*theta;

while ((sqrt(update'*update) > theta) & (iter < Max_iter))
   iter = iter + 1;
   
   %Find all incorrectly classified samples, Yk
   Yk		= find(a'*train_patterns.*(2*train_targets-1) < 0);
   update	= eta * sum(y(:,Yk)')';
   
   % a <- a + eta*sum(Yk)
   a = a + update;
   
   % Save updates 
   updates(iter) = sum(abs(update));
end

if (iter == theta),
   disp(['Maximum iteration (' num2str(theta) ') reached']);
end

%Classify test patterns
test_targets = a'*[test_patterns; ones(1, size(test_patterns,2))] > 0;