function [test_targets, errors] = Components_without_DF(train_patterns, train_targets, test_patterns, Classifiers)

% Classify points using component classifiers without discriminant functions
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	Classifiers		- Classification algorithms. The format is:
%							[('<Algorithm 1>', <Algorithm 1 parameters>), ...]
%
% Outputs
%	test_targets	- Predicted targets
%   errors          - Errors throughout the iterations
%

%Read the algorithms
open_bracket	= findstr(Classifiers,'(');
close_bracket	= findstr(Classifiers,')');
if length(open_bracket) ~= length(close_bracket),
    error('Input vector contains an error!')
end
Nalgorithms		= length(open_bracket);
for i = 1:Nalgorithms,
    line	= Classifiers(open_bracket(i)+1:close_bracket(i)-1);
    comma = findstr(line,',');
    if isempty(comma),
        algorithms(i).name      = line(2:end-1);
        algorithms(i).params    = [];
    else
        algorithms(i).name      = line(2:comma-2);
        algorithms(i).params	= str2num(line(comma+1:end));
    end
end

%Train the weak classifiers
disp('Training weak classifiers')
[Ndim, M]   		= size(train_patterns);
p				    = zeros(Nalgorithms, M);
weak_test_targets   = zeros(Nalgorithms, size(test_patterns,2));
weak_train_targets  = zeros(Nalgorithms, M);
for i = 1:Nalgorithms,
    weak_targets            = feval(algorithms(i).name, train_patterns, train_targets, [train_patterns, test_patterns], algorithms(i).params);
    weak_train_targets(i,:) = weak_targets(1:M);
    weak_test_targets(i,:)  = weak_targets(M+1:end);
end

p   = exp(weak_train_targets)./(1+exp(1));   %Use the softmax transformation of the data. We only have {0,1} classes, so the transformation is simple

%Init gating components
Ndim			= Ndim + 1;
x				= [train_patterns; ones(1,M)];
y				= train_targets;
alpha			= randn(Ndim, Nalgorithms);
alpha			= sqrtm(cov(x',1)+randn(Ndim))*alpha + mean(x')'*ones(1,Nalgorithms);
w   			= exp(alpha'*x)./(ones(Nalgorithms,1)*sum(exp(alpha'*x)));

%Learn the gating parameters
disp('Finding gating parameters')
errors			= 1e3;
derror          = 1;
iter            = 1;

while ((errors(iter) > 1/M) & (derror > 0)),    
    iter = iter + 1;
    
    %Perform gradient descent on the gating parameters
    h               = w.*p./(ones(Nalgorithms,1)*sum(w.*p));
    dalpha          = (x*(h - w)');
    alpha			= alpha + dalpha;
    
    w				= exp(alpha'*x)./(ones(Nalgorithms,1)*sum(exp(alpha'*x)));
    Y				= sum(w.*p);
    errors(iter)    = sum(y ~= (Y>.5))/M;
    
    derror          = errors(iter) - errors(iter-1);
    disp(['Error is ' num2str(errors(iter))]) 
end

%Classify test patterns
test_patterns   = [test_patterns; ones(1, size(test_patterns,2))];
u			    = exp(alpha'*test_patterns)./(ones(Nalgorithms,1)*sum(exp(alpha'*test_patterns)));
test_targets    = sum(u.*(exp(weak_test_targets)/(1+exp(1))));

if (length(unique(train_targets)) == 2)
    test_targets = test_targets > 0.5;
end