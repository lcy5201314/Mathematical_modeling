function [test_targets, Wh, Wo, J] = Optimal_Brain_Surgeon(train_patterns, train_targets, test_patterns, params)

% Classify using a backpropagation network with a batch learning algorithm and remove excess units
% using the optimal brain surgeon algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	params          - Initial number of hidden units, Convergence criterion
%
% Outputs
%	test_targets	- Predicted targets
%   Wh              - Hidden unit weights
%   Wo              - Output unit weights
%   J               - Errors throughout the training

[Nh, Theta] = process_params(params);
[Ni, M]		= size(train_patterns);
Uc          = unique(train_targets);

%Train a reasonably large network to minimum error
disp(['Training a neural network with ' num2str(Nh) ' hidden units']);
[test_targets, Wh, Wo] = Backpropagation_Batch(train_patterns, train_targets, test_patterns, [Nh, Theta, 0.1]);

if (length(Uc) == 2)
    train_targets = (train_targets>0)*2-1;
end

J             = 100;
gradJ         = 0;
iter          = 1;

%Cascade all the weights
W       = [Wo(:); Wh(:)];

%Define initial H's
alpha   = 1e-8;
invH    = (alpha^-1)*eye(length(W));         
Lq      = [1 0];
pruned  = 1:length(W);

disp('Pruning excess units');
while ((gradJ < Theta) & (length(unique(Lq)) > 1)),

    %Compute inv(H) by equation 70, DHS chapter 6
    for i = 1:M,
        xm = train_patterns(:,i);
        tk = train_targets(i);
        
        %Forward propagate the input:
        %First to the hidden units
        gh				= Wh*[xm; 1];
        [y, dfh]		= activation(gh);
        
        %Now to the output unit
        go				= Wo*[y; 1];
        [zk, dfo]	    = activation(go);
  
        Ni              = size(xm,1)+1;
        Xv              = dfo*[y; 1];
        Xu              = dfo*(dfh*ones(1,Ni)).*Wh.*(y*ones(1,Ni));

        Xm              = [Xv; Xu(:)];

        invH            = invH + (invH*Xm*Xm'*invH)/(M + Xm'*invH*Xm);
    end
    
    if ~isfinite(sum(sum(invH))),
        break
    end
    
    %q* <- argmin(w_q^2/(2*inv(H)_qq))
    Lq          = W.^2./(2*diag(invH));
    [m, q_star] = min(Lq(pruned)); %We don't want to prune the same weight twice
    q_star      = pruned(q_star);
    e_q_star    = zeros(length(W), 1);
    e_q_star(q_star) = 1;
    
    pruned = pruned(find(pruned ~= q_star));
    
    %w  <- w - w*_q/inv(H)_q*q**inv(H)*e_q*
    W           = W - W(q_star)/invH(q_star, q_star)*invH*e_q_star;
    
    Wo          = W(1:size(Wo,2))';
    Wh          = reshape(W(size(Wo,2)+1:end),size(Wh));

    iter = iter + 1;
    %Calculate total error
    J(iter) = 0;
    for i = 1:M,
        J(iter) = J(iter) + (train_targets(i) - activation(Wo*[activation(Wh*[train_patterns(:,i); 1]); 1])).^2;
    end
    J(iter) = J(iter)/M; 
    gradJ = J(iter) - J(iter-1);

    disp(['Removed weight number ' num2str(q_star) '. Gradient jump was ' num2str(gradJ)]);

end

%Compute the test targets
test_targets = zeros(1, size(test_patterns,2));
for i = 1:size(test_patterns,2),
    Xm              = test_patterns(:,i);
    test_targets(i) = activation(Wo*[activation(Wh*[Xm; 1]); 1]);
end

%If there are two classes, collapse the output values
if (length(Uc) == 2)
    test_targets = test_targets > 0;
end


function [f, df] = activation(x)
%The activation function for a neural network
a = 1.716;
b = 2/3;
f	= a*tanh(b*x);
df	= a*b*sech(b*x).^2;