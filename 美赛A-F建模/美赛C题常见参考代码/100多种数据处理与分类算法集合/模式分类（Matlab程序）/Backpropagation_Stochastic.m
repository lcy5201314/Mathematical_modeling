function [test_targets, Wh, Wo, J] = Backpropagation_Stochastic(train_patterns, train_targets, test_patterns, params)

% Classify using a backpropagation network with stochastic learning algorithm
% Inputs:
% 	training_patterns   - Train patterns
%	training_targets	- Train targets
%   test_patterns       - Test  patterns
%	params              - Number of hidden units, Convergence criterion, Convergence rate
%
% Outputs
%	test_targets        - Predicted targets
%   Wh                  - Hidden unit weights
%   Wo                  - Output unit weights
%   J                   - Error throughout the training

[Nh, Theta, eta] = process_params(params);
iter	         = 1;

[Ni, M]          = size(train_patterns);
No		         = 1;

Uc               = length(unique(train_targets));
%If there are only two classes, remap to {-1,1}
if (Uc == 2)
    train_targets    = (train_targets>0)*2-1;
end

%Initialize the net: In this implementation there is only one output unit, so there
%will be a weight vector from the hidden units to the output units, and a weight matrix
%from the input units to the hidden units.
%The matrices are defined with one more weight so that there will be a bias
w0		= max(abs(std(train_patterns')'));
Wh		= rand(Nh, Ni+1).*w0*2-w0; %Hidden weights
Wo		= rand(No, Nh+1).*w0*2-w0; %Output weights

Wo    = Wo/mean(std(Wo'))*(Nh+1)^(-0.5);
Wh    = Wh/mean(std(Wh'))*(Ni+1)^(-0.5);

rate	= 10*Theta;
J(1)    = 1e3;

while (rate > Theta),
    %Randomally choose an example
    i	= randperm(M);
    m	= i(1);
    Xm = train_patterns(:,m);
    tk = train_targets(m);
    
    %Forward propagate the input:
    %First to the hidden units
    gh				= Wh*[Xm; 1];
    [y, dfh]		= activation(gh);
    %Now to the output unit
    go				= Wo*[y; 1];
    [zk, dfo]	= activation(go);
    
    %Now, evaluate delta_k at the output: delta_k = (tk-zk)*f'(net)
    delta_k		= (tk - zk).*dfo;
    
    %...and delta_j: delta_j = f'(net)*w_j*delta_k
    delta_j		= dfh'.*Wo(1:end-1).*delta_k;
    
    %w_kj <- w_kj + eta*delta_k*y_j
    Wo				= Wo + eta*delta_k*[y;1]';
    
    %w_ji <- w_ji + eta*delta_j*[Xm;1]
    Wh				= Wh + eta*delta_j'*[Xm;1]';
    
    iter 			= iter + 1;

    %Calculate total error
    J(iter)    = 0;
    for i = 1:M,
        J(iter) = J(iter) + (train_targets(i) - activation(Wo*[activation(Wh*[train_patterns(:,i); 1]); 1])).^2;
    end
    J(iter) = J(iter)/M; 
    rate  = abs(J(iter) - J(iter-1))/J(iter-1)*100;
    
    if (iter/100 == floor(iter/100)),
        disp(['Iteration ' num2str(iter) ': Total error is ' num2str(J(iter))])
    end
    
end

disp(['Backpropagation converged after ' num2str(iter) ' iterations.'])

%Classify the test patterns
test_targets = zeros(1, size(test_patterns,2));
for i = 1:size(test_patterns,2),
    test_targets(i) = activation(Wo*[activation(Wh*[test_patterns(:,i); 1]); 1]);
end

if (Uc == 2)
    test_targets  = test_targets >0;
end



function [f, df] = activation(x)

a = 1.716;
b = 2/3;
f	= a*tanh(b*x);
df	= a*b*sech(b*x).^2;