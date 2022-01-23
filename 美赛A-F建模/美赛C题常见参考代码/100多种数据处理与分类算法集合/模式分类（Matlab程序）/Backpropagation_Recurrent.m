function [test_targets, W, J] = Backpropagation_Recurrent(train_patterns, train_targets, test_patterns, params)

% Classify using a backpropagation recurrent network with a batch learning algorithm
% Inputs:
% 	training_patterns   - Train patterns
%	training_targets	- Train targets
%   test_patterns       - Test  patterns
%	params              - Number of hidden units, Convergence criterion, Convergence rate
%
% Outputs
%	test_targets        - Predicted targets
%   W                   - Unit weights
%   J                   - Error throughout the training

[Nh, Theta, eta]    = process_params(params);
[Ni, M] 	        = size(train_patterns);
eta		            = eta/M;
iter		        = 1;
NiterDisp           = 1;
maxIter 	        = 100;
Uc                  = length(unique(train_targets));

No		            = 1; %Number of output units
Ni		    	    = Ni + 1;
train_patterns      = [train_patterns; ones(1,M)];

%Init the network weights
w0		= max(abs(std(train_patterns')'));
W		= rand(Nh+No, Ni+No+Nh).*w0*2-w0; %The first Nh units are hidden, and the last No are the output
W       = W/mean(std(W'))*(Nh+No)^(-0.5);

%If there are only two classes, remap to {-1,1}
if (Uc == 2)
    train_targets    = (train_targets>0)*2-1;
end

rate	= 10*Theta;
J(1)    = 1e3;

while (rate > Theta),
    %Randomally choose an example
    i	= randperm(M);
    m	= i(1);
    Xm = train_patterns(:,m);
    tk = train_targets(m);

    dOut    = 1;   %Fraction of change in output at time t
    dOutOld = 2;
    y_k     = zeros(Nh+No, 1);
    P_k     = zeros(Nh+No, Nh+No, Ni+Nh+No);
    U       = zeros(Nh+No, Nh+No, Ni+Nh+No);
    C       = [zeros(1, Nh) ones(1, No)];
    deltaW  = zeros(size(W));
    iter1    = 0;
    
    %Now, introduce the input to the net, and continue to feed it until it stabilizes
    while ((abs(dOut-dOutOld)/dOutOld > 1e-3) & (iter1 < maxIter)),
        iter1 = iter1 + 1;
        
        %Compute network output
        z_k                 = [Xm; y_k];
        net_k               = W*z_k;
        [y_k_tplus1, dy_k]  = activation(net_k);
        
        %Compute the error
        e_k                 = (tk - C*y_k_tplus1)';
        
        dOutOld             = dOut;
        dOut                = sum(abs(e_k./tk));
        
        %Build U and Phi
        for i = 1:Nh+No,
            U(i,i,:) = z_k';
        end
        Phi = eye(No+Nh) .* (dy_k*ones(1,No+Nh));
        
        %Compute the weight update
        for i = 1:Nh+No,
            deltaW(i,:)     = deltaW(i,:) + eta * C * squeeze(P_k(i,:,:)) * e_k;
        end
                
        %Update P_k
        for i = 1:Nh+No,
            P_k(i,:,:)      = Phi * (W(:,Ni+1:end) * squeeze(P_k(i,:,:)) + squeeze(U(i,:,:)));
        end

        y_k                 = y_k_tplus1;

    end

    %Update the weights
    W                   = W + deltaW;

    iter 			= iter + 1;

    %Measure the error
    J(iter) = 0;
    for i = 1:M,
        Xm = train_patterns(:,i);
        tk = train_targets(i);

        dOut    = 1;   %Fraction of change in output at time t
        dOutOld = 2;
        y_k     = zeros(Nh+No, 1);
        iter1    = 0;
        
       %Now, introduce the input to the net, and continue to feed it until it stabilizes
        while ((abs(dOut-dOutOld)/dOutOld > 1e-3) & (iter1 < maxIter)),
            iter1 = iter1 + 1;
            %Compute network output
            z_k     = [Xm; y_k];
            net_k   = W*z_k;
            y_k     = activation(net_k);
            e_k     = (tk - C*y_k)';
            dOutOld = dOut;
            dOut    = sum(abs((tk-e_k)./tk));
        end
        
        J(iter)   = J(iter) + (tk - e_k).^2;        
    end
    J(iter) = J(iter)/M; 
    rate  = abs(J(iter) - J(iter-1))/J(iter-1)*100;
               
    if (iter/NiterDisp == floor(iter/NiterDisp)),
        disp(['Iteration ' num2str(iter) ': Error is ' num2str(J(iter))])
    end
    
end

disp(['Backpropagation converged after ' num2str(iter) ' iterations.'])

%Classify the test patterns
test_targets = zeros(1, size(test_patterns,2));

for i = 1:size(test_patterns,2),
    Xm = [test_patterns(:,i); 1];
        
    dOut    = 1;   %Fraction of change in output at time t
    dOutOld = 2;
    y_k     = zeros(Nh+No, 1);
    iter1    = 0;
    
    %Now, introduce the input to the net, and continue to feed it until it stabilizes
    while ((abs(dOut-dOutOld)/dOutOld > 1e-3) & (iter1 < maxIter)),
        iter1 = iter1 + 1;
        
        %Compute network output
        z_k     = [Xm; y_k];
        net_k   = W*z_k;
        y_k     = activation(net_k);
        e_k     = (tk - C*y_k);
        dOutOld = dOut;
        dOut    = sum(abs((tk-e_k)./tk));
    end

    test_targets(i) = y_k(end);
end

if (Uc == 2)
    test_targets  = test_targets >0;
end



function [f, df] = activation(x)

a = 1.716;
b = 2/3;
f	= a*tanh(b*x);
df	= a*b*sech(b*x).^2;