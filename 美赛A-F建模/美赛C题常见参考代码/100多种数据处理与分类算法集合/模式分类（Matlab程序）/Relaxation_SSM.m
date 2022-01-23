function [test_targets, a] = Relaxation_SSM(train_patterns, train_targets, test_patterns, params)

% Classify using the single-sample relaxation with margin algorithm
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

[c, n]		    = size(train_patterns);
[Max_iter, b, eta]	 = process_params(params);

y               = [train_patterns ; ones(1,n)];
train_zero      = find(train_targets == 0);

%Preprocessing
processed_patterns = y;
processed_patterns(:,train_zero) = -processed_patterns(:,train_zero);

%Initial weights
a               = sum(processed_patterns')';
iter  	        = 0;
k				= 0;

while ((sum(a'*processed_patterns < b)>0) & (iter < Max_iter))
    iter = iter + 1;
    
    %k <- (k+1) mod n
    k = mod(k+1,n);
    if (k == 0), 
        k = n;
    end
    
    if (a'*processed_patterns(:,k) <= b),
        % a <- a + eta*sum((b-w'*Yk)/||Yk||*Yk)
        grad			= (b-a'*y(:,k))./sum(y(:,k).^2);
        update		= grad.*y(:,k);
        a 	= a + eta * update;
    end
    
end

if (iter == Max_iter),
    disp(['Maximum iteration (' num2str(Max_iter) ') reached']);
end

%Classify test patterns
test_targets = a'*[test_patterns; ones(1, size(test_patterns,2))] > 0;