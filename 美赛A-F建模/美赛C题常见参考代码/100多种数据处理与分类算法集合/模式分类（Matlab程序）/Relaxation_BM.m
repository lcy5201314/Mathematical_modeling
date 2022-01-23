function [test_targets, a] = Relaxation_BM(train_patterns, train_targets, test_patterns, params)

% Classify using the batch relaxation with margin algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	param		    - [Max iter, Margin, Convergence rate]
%
% Outputs
%	test_targets	- Predicted targets
%   a               - Classifier weights
%
% NOTE: Works for only two classes.

[c, n]				= size(train_patterns);
[Max_iter, b, eta]	= process_params(params);

y               = [train_patterns ; ones(1,n)];
train_zero      = find(train_targets == 0);

%Preprocessing
processed_patterns = y;
processed_patterns(:,train_zero) = -processed_patterns(:,train_zero);

%Initial weights
a               = sum(processed_patterns')';
iter  	        = 0;
Yk				= [1];

while (~isempty(Yk) & (iter < Max_iter))
   iter = iter + 1;
   
   %If a'y_j <= b then append y_j to Yk
   Yk = [];
   for k = 1:n,
   	if (a'*processed_patterns(:,k) <= b),
         Yk = [Yk k];
      end
   end
   
   if isempty(Yk),
      break
   end
   
   % a <- a + eta*sum((b-w'*Yk)/||Yk||*Yk)
   grad			= (b-a'*y(:,Yk))./sum(y(:,Yk).^2);
   update		= sum(((ones(c+1,1)*grad).*y(:,Yk))')';
   a            = a + eta * update;
end

if (iter == Max_iter),
   disp(['Maximum iteration (' num2str(Max_iter) ') reached']);
end

%Classify test patterns
test_targets = a'*[test_patterns; ones(1, size(test_patterns,2))] > 0;