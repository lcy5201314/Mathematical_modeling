function [test_targets, a] = Perceptron(train_patterns, train_targets, test_patterns, alg_param)

% Classify using the Perceptron algorithm (Fixed increment single-sample perceptron)
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	alg_param	    - Either: Number of iterations, weights vector or [weights, number of iterations]
%
% Outputs
%	test_targets	- Predicted targets
%   a               - Perceptron weights
%
% NOTE: Works for only two classes.

[c, r]		   = size(train_patterns);

%Weighted Perceptron or not?
switch length(alg_param),
case r + 1,
    %Ada boost form
    p           = alg_param(1:end-1);
    max_iter    = alg_param(end);
case {r,0},
    %No parameter given
    p           = ones(1,r);
    max_iter    = 5000;
otherwise
    %Number of iterations given
    max_iter    = alg_param;
    p           = ones(1,r);
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
   indice = 1 + floor(rand(1)*n);
   if (a' *  y(:,indice) <= 0)
      a = a + p(indice)* y(:,indice);
   end
end

if (iter == max_iter)&(length(alg_param)~= r + 1),
   disp(['Maximum iteration (' num2str(max_iter) ') reached']);
end

%Classify test patterns
test_targets = a'*[test_patterns; ones(1, size(test_patterns,2))] > 0;