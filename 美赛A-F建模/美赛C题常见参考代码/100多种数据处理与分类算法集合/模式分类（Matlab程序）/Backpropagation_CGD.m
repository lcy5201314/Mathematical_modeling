function [test_targets, Wh, Wo, errors] = Backpropagation_CGD(train_patterns, train_targets, test_patterns, params)

% Classify using a backpropagation network with a batch learning algorithm and conjugate gradient descent
% Inputs:
% 	training_patterns   - Train patterns
%	training_targets	- Train targets
%   test_patterns       - Test  patterns
%	params              - Number of hidden units, Convergence criterion
%
% Outputs
%	test_targets        - Predicted targets
%   Wh                  - Hidden unit weights
%   Wo                  - Output unit weights
%   J                   - Error throughout the training

[Nh, Theta] = process_params(params);
iter	    = 0;
IterDisp    = 5;

[Ni, M]     = size(train_patterns);
No          = 1;
Uc          = length(unique(train_targets));

%If there are only two classes, remap to {-1,1}
if (Uc == 2)
    train_targets    = (train_targets>0)*2-1;
end

%Initialize the net: In this implementation there is only one output unit, so there
%will be a weight vector from the hidden units to the output units, and a weight matrix
%from the input units to the hidden units.
%The matrices are defined with one more weight so that there will be a bias
w0			= max(abs(std(train_patterns')'));
Wh			= rand(Nh, Ni+1).*w0*2-w0; %Hidden weights
Wo			= rand(No,  Nh+1).*w0*2-w0; %Output weights
Wo          = Wo/mean(std(Wo'))*(Nh+1)^(-0.5);
Wh          = Wh/mean(std(Wh'))*(Ni+1)^(-0.5);

%Iteration zero: Compute gradJ and from it the CGD matrices
deltaWo	= 0;
deltaWh	= 0;

for m = 1:M,
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
    
    %delta_w_kj <- w_kj + eta*delta_k*y_j
    deltaWo		= deltaWo + delta_k*[y;1]';
    
    %delta_w_ji <- w_ji + eta*delta_j*[Xm;1]
    deltaWh		= deltaWh + delta_j'*[Xm;1]';
    
end
S       = [deltaWh(:); deltaWo(:)];
R       = [deltaWh(:); deltaWo(:)];
oldR    = zeros(size(R));
normR0  = sqrt(sum(R.^2));

%Now, for each iteration
while 1,
    %Use a line search to find eta that minimizes J(eta)
    Eta = logspace(-5, 1, 50);
    for i = 1:length(Eta),
        J(i)  = find_error(Eta(i), Wo, reshape(S(Nh*(Ni+1)+1:end), No, Nh+1), Wh, reshape(S(1:Nh*(Ni+1)), Nh, Ni+1), train_patterns, train_targets);
    end
    [m, i]  = min(J);
    eta     = Eta(i);
    
    %Update iteration number
    iter 			= iter + 1;
    if (iter/IterDisp == floor(iter/IterDisp)),
        disp(['Iteration ' num2str(iter) ': Total error is ' num2str(J(i)) '. Step size is: ' num2str(eta)])
    end
    
    errors(iter) = sqrt(sum(R.^2));
    
    %Test if the residual has decreased enough
    if (sqrt(sum(R.^2)) < Theta*normR0),
        break
    end
    
    %Update the weight matrices
    Wo  = Wo + eta * reshape(S(Nh*(Ni+1)+1:end), No, Nh+1);
    Wh  = Wh + eta * reshape(S(1:Nh*(Ni+1)), Nh, Ni+1);
    
    %Compute the new gradient matrices using backpropagation
    deltaWo	= 0;
    deltaWh	= 0;
    
    for m = 1:M,
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
        
        %delta_w_kj <- w_kj + eta*delta_k*y_j
        deltaWo		= deltaWo + delta_k*[y;1]';
        
        %delta_w_ji <- w_ji + eta*delta_j*[Xm;1]
        deltaWh		= deltaWh + delta_j'*[Xm;1]';
        
    end
    
    %Set R
    R = [deltaWh(:); deltaWo(:)];
    
    %Use the Polak-Ribiere method to calculate beta
    if ((oldR'*oldR) ~= 0),
        beta  = max(0, R'*(R - oldR)/(oldR'*oldR));
    else
        beta  = 0;
    end
    
    %Update the direction vector
    S       = R + beta * S;
    
    %Update the old vectors
    oldR    = R;
    
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

function J = find_error(eta, Wo, deltaWo, Wh, deltaWh, patterns, targets)
%Find the error given the net parameters
M    = size(patterns,2);
J    = 0;
for i = 1:M,
    J = J + ((targets(i) - activation((Wo+eta*deltaWo)*[activation((Wh+eta*deltaWh)*[patterns(:,i); 1]); 1])).^2);
end
J    = J/M;



function [f, df] = activation(x)

a = 1.716;
b = 2/3;
f	= a*tanh(b*x);
df	= a*b*sech(b*x).^2;