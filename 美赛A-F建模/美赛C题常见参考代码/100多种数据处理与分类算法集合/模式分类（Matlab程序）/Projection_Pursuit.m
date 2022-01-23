function [test_targets, V, Wo] = Projection_Pursuit(train_patterns, train_targets, test_patterns, Ncomponents)

% Classify using projection pursuit regression
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	Ncomponents     - Number of components to project on
%
% Outputs
%	test_targets	- Predicted targets
%   V               - Components weights
%   Wo              - Output unit weights

[Ni, M] = size(train_patterns);

iter			= 1;
NiterDisp	    = 1;
maxIter 		= 50;

train_patterns= [train_patterns; ones(1,M)];
Ni		  = Ni + 1;
No		  = 1; %Number of output units

V		  = rand(Ni, Ncomponents);
Wo		  = rand(No, Ncomponents+1);
gradJ	  = 1;
J         = 1000;

%Find the regression parameters
while ((gradJ > 1e-2) & (iter < maxIter)),
   iter = iter + 1;
   
   %Optimize for the components
   J1 = inline('sum((t - Wo*[tanh(V''*x); ones(No, M)]).^2)','V','t','Wo','x','No','M');
	V  = fminunc(J1, V, [], train_targets, Wo, train_patterns, No, M);       
   
   %Optimize the weights
   J2 = inline('sum((t - Wo*[tanh(V''*x); ones(No, M)]).^2)','Wo','t','V','x','No','M');
	Wo = fminunc(J2, Wo, [], train_targets, V, train_patterns, No, M);       
   
   %Evaluate the error
   J(iter) = feval(J1, V, train_targets, Wo, train_patterns, No, M);
   gradJ   = abs(J(iter) - J(iter-1));
end

if (iter == maxIter),
   disp('Optimization terminated after reaching the maximum iteration.')
else
   disp(['Converged after ' num2str(iter) ' iterations.'])
end

%Classify test patterns
N            = size(test_patterns,2);
test_targets = Wo*[tanh(V'*[test_patterns; ones(1, N)]); ones(No, N)];

%If there are only two classes, collapse them
if (length(unique(train_targets)) == 2)
    test_targets = test_targets > 0.5;
end
