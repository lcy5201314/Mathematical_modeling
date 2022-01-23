function [test_targets, w_percept, b] = Ho_Kashyap(train_patterns, train_targets, test_patterns, params)

% Classify using the using the Ho-Kashyap algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	params		    - [Type(Basic/Modified), Maximum iteration, Convergence criterion, Convergence rate]
%
% Outputs
%	test_targets	- Predicted targets
%   w_percept       - Classifier weights
%   b               - Margin

[c, n]		   = size(train_patterns);

[type, Max_iter, b_min, eta] = process_params(params);

train_patterns  = [train_patterns ; ones(1,n)];
train_zero      = find(train_targets == 0);

%Preprocessing (Needed so that b>0 for all patterns)
processed_patterns = train_patterns;
processed_patterns(:,train_zero) = -processed_patterns(:,train_zero);
processed_targets  = 2*train_targets-1;

b                  = ones(1,n);
Y                  = processed_patterns;
a                  = pinv(Y')*b';
k	               = 0;
e    	           = 1e3;
found              = 0;

while ((sum(abs(e) > b_min)>0) & (k < Max_iter) &(~found))

    %k <- (k+1) mod n
    k = k+1;

    %e <- Ya - b
    e       = (Y' * a)' - b;
   
    %e_plus <- 1/2(e+abs(e))
    e_plus  = 0.5*(e + abs(e));
   
    %b <- b + 2*eta*e_plus
    b       = b + 2*eta*e_plus;
    
    if strcmp(type,'Basic'),
        %a <- pinv(Y)*b;   
        a = pinv(Y')*b';
    else
        %a <- a + eta*pinv(Y)*|e_plus|;   
        a = a + eta*pinv(Y')*abs(e_plus)';
    end        
    
end

if (k == Max_iter),
   disp(['No solution found']);
else
   disp(['Did ' num2str(k) ' iterations'])
end

test_targets = a' * [test_patterns; ones(1, size(test_patterns,2))];

if (length(unique(train_targets)) == 2)
    test_targets = test_targets > 0;
end