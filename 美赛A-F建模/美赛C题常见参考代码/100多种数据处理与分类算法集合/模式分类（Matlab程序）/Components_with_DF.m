function [test_targets, errors] = Components_with_DF(train_patterns, train_targets, test_patterns, Ncomponents)

% Classify points using component classifiers with discriminant functions
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	Ncomponents		- Number of component classifiers
%
% Outputs
%	test_targets	- Predicted targets
%   errors          - Errors throughout the iterations
%
% This implementation works with logistic component classifiers and a softmax gating function
% The parameters of the components are learned using Newton descent, and the parameters
% of the gating system using gradient descent

[Ndim, M] 	    = size(train_patterns);
Ndim			= Ndim + 1;
x				= [train_patterns; ones(1,M)];
y				= train_targets;
theta			= zeros(Ndim, Ncomponents)+eps;
alpha			= randn(Ndim, Ncomponents);
alpha			= sqrtm(cov(x',1)+randn(Ndim))*alpha + mean(x')'*ones(1,Ncomponents);

derror          = 1000;
errors			= 1;
iter            = 1;

while ((errors(iter) > 1/M) & (derror > 0)),
    
    %Perform gradient descent on the component classifiers
    w			= exp(alpha'*x)./(ones(Ncomponents,1)*sum(exp(alpha'*x)));
    for i = 1:Ncomponents,
        p					= exp(theta(:,i)'*x)./(1+exp(theta(:,i)'*x));
        W					= diag(p.*(1-p));
        delta_theta_i	= inv(x*W*x')*x*(y.*w(i,:) - p)';
        if ~isfinite(sum(delta_theta_i)),
            delta_theta_i = 0;
        end
        theta(:,i)		= theta(:,i) + delta_theta_i;
    end
    
    %Perform gradient descent on the gating parameters
    p				= zeros(Ncomponents, M);
    for i = 1:Ncomponents,
        p(i,:)			= exp(theta(:,i)'*x)./(1+exp(theta(:,i)'*x));
    end
    h               = w.*p./(ones(Ncomponents,1)*sum(w.*p));
    dalpha          = (x*(h - w)');
    alpha			= alpha + dalpha;
    
    iter            = iter + 1;
    w				= exp(alpha'*x)./(ones(Ncomponents,1)*sum(exp(alpha'*x)));
    Y				= sum(w.*p);
    errors(iter)	= sum(y ~= (Y>.5))/M;
    
    derror          = errors(iter) - errors(iter-1);
    disp(['Error is ' num2str(errors(iter))]) 
end

%Classify test patterns
test_patterns   = [test_patterns; ones(1,size(test_patterns,2))];
y				= exp(theta'*test_patterns)./(ones(Ncomponents,size(test_patterns,2)) + exp(theta'*test_patterns));
u				= exp(alpha'*test_patterns)./(ones(Ncomponents,1)*sum(exp(alpha'*test_patterns)));
test_targets    = sum(y.*u);

if (length(unique(train_targets)) == 2)
    test_targets = test_targets > 0.5;
end