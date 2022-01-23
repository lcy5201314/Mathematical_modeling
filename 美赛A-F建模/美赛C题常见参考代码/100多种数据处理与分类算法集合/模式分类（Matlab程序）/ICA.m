function [patterns, targets, W, Aw, means] = ICA(patterns, targets, params)

%Reshape the data points using the independent component analysis algorithm
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	params			- [Output dimension, Learning rate]
%
%Outputs
%	patterns		- New patterns
%	targets			- New targets
%	W				- Reshape martix
%   means           - The means vector of the patterns              

[r,c]			= size(patterns);
[dimension, eta] = process_params(params);

if (r < dimension),
    error('Output dimension cannot be larger than the input dimension')
end

%Whiten the data to zero mean and unit covariance
means       = mean(patterns')';
patterns    = patterns - means*ones(1,c);
[v, d]	    = eig(cov(patterns',1));
Aw			= v*inv(sqrtm(d));
patterns    = Aw'*patterns;

%Move data to the range of [-1,1]
minp     = min(patterns')';
maxp     = max(patterns')';
patterns = (patterns - minp*ones(1,c))./((maxp-minp)*ones(1,c));
patterns = patterns*2-1;

%Find the weight matrix
W			= randn(r);
iter        = 1;
while (iter < 1000),
    iter    = iter + 1;
    y		= W*patterns;
    phi	    = activation(y);
    dW		= (eye(r) - 1/c*phi*y')*W;
    
    %Break if algorithm diverges
    if (max(max(dW)) > 1e3),
        disp(['Algorithm diverged after ' num2str(i) ' iterations'])
        break
    end
    
    W		= W + eta*dW;   
    
    update  = max(max(abs(dW)));
    
    %If the algorithm converged, exit
    if (update < eta),
        disp(['Algorithm converged after ' num2str(iter) ' iterations'])
        break
    else
        if (iter / 10 == floor(iter/10))
            disp(['Iteration ' num2str(iter) ': Maximum update is ' num2str(update)])
        end
    end
end

%Take only the most influential outputs
power		= sum(abs(W)');
[m, in]	    = sort(power);
W			= W(in(r-dimension+1:r),:);

%Calculate new patterns
patterns = W*patterns;

W        = W*Aw;

%End ICA

function phi = activation(y)
%Activation function for ICA
phi = y.^3;
%phi = tahn(y);