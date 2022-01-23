function [test_targets, w_pocket] = Pocket(train_patterns, train_targets, test_patterns, alg_param)

% Classify using the pocket algorithm (an improvement on the perceptron)
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	alg_param	    - Either: Number of iterations, weights vector or [weights, number of iterations]
%
% Outputs
%	test_targets	- Predicted targets
%   w               - Pocket weights
%
% NOTE: Works for only two classes.

[c, r] = size(train_patterns);

%Weighted Pocket or not?
switch length(alg_param),
case r + 1,
    %Ada boost form
    p           = alg_param(1:end-1);
    max_iter    = alg_param(end);
case {r, 0},
    %No parameter given
    p           = ones(1,r);
    max_iter    = 500;
otherwise
    %Number of iterations given
    max_iter    = alg_param;
    p           = ones(1,r);
end

train_patterns = [train_patterns ; ones(1,r)];
train_one      = find(train_targets == 1);
train_zero     = find(train_targets == 0);

%Preprocessing
processed_patterns = train_patterns;
processed_patterns(:,train_zero) = -processed_patterns(:,train_zero);

%Initial weights
w_percept   = sum(processed_patterns')';
%w_percept   = train_patterns .* (ones(c+1,1) * (2*(train_targets-0.5)));
%w_percept	= rand(c+1,1);
w_pocket	= rand(c+1,1);

correct_classified = 0;
n						 = length(train_targets);
iter					 = 0;

while ((longest_run(w_percept, processed_patterns) < n) & (iter < max_iter))
   iter = iter + 1;
   %Every 10 points, do the pocket switchover
   for i = 1:10,
      indice = 1 + floor(rand(1)*n);
      if (w_percept' * processed_patterns(:,indice) <= 0)
         w_percept = w_percept + p(indice) * processed_patterns(:,indice);
      end
   end
   %Find if it is neccessary to change weights:
   if (longest_run(w_percept, processed_patterns) > longest_run(w_pocket, processed_patterns)),
      w_pocket = w_percept;
   end
end

if (iter == max_iter)&(length(alg_param)~= r + 1),
   disp(['Maximum iteration (' num2str(max_iter) ') reached']);
end

%Classify test patterns
test_targets = w_pocket'*[test_patterns; ones(1, size(test_patterns,2))] > 0;

%END

function L = longest_run(weights, patterns)

%Find the length of the longest run of correctly classified random points
n           = length(patterns);
indices     = randperm(n);
L           = 0;
correct     = 1;

for i = 1:n,
   if (weights' * patterns(:,indices(i)) <= 0)	%Find if it is correctly classified
      break
   end
   L = i;
end
