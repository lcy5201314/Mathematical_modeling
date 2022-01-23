function [test_targets, a] = Perceptron_VIM(train_patterns, train_targets, test_patterns, params)

% Classify using the variable incerement Perceptron with margin algorithm 
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	param	    	- [Num iter, Margin, Convergence rate]
%
% Outputs
%	test_targets	- Predicted targets
%   a               - Perceptron weights
%
% NOTE: Works for only two classes.

[c, n]	 	       = size(train_patterns);
[Max_iter, b, eta] = process_params(params);

train_patterns  = [train_patterns ; ones(1,n)];
train_zero      = find(train_targets == 0);

%Preprocessing
y               = train_patterns;
y(:,train_zero) = -y(:,train_zero);

%Initial weights
a               = sum(y')';

iter			= 0;
k				= 0;

while ((sum(a'*y <= b)>0) & (iter < Max_iter))
   iter = iter + 1;
   
   %k <- (k+1) mod n
   k = mod(k+1,n);
   if (k == 0), 
      k = n;
   end
   
   if (a'*y(:,k) <= b),
	   % a <- a + eta*sum(Yk)
      a = a + eta * y(:,k);
   end
   
end

if (iter == Max_iter),
   disp(['Maximum iteration (' num2str(Max_iter) ') reached']);
end

%Classify test patterns
test_targets = a'*[test_patterns; ones(1, size(test_patterns,2))] > 0;