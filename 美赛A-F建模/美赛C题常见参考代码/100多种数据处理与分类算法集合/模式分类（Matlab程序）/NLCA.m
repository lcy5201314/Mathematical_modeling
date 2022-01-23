function [new_patterns, targets] = NLPCA(patterns, targets, params)

% Reshape the data using the non-linear PCA algorithm
% Inputs:
% 	patterns - Train patterns
%	targets	 - Train targets
%	params   - [Number of output dimensions, Number of hidden units]
%
% Outputs
%	patterns - New patterns
%	targets  - Targets

[Ni, M]    = size(patterns);
[Dout, Nh] = process_params(params);

Theta   = 0.05;
eta     = 0.1/M;
iter	= 0;
IterDisp= 1;

patterns= [patterns; ones(1,M)];
Ni		= Ni + 1;

if (Dout > Ni),
   error('Error: The output dimension should be smaller than the input dimension')
end

%Initialize the net: 
%This net has two non-linear hidden layers separated by a linear layer
%The output layer is linear
w0		= max(abs(std(patterns')'));
Wnl1	= rand(Nh, Ni).*w0*2-w0; %Hidden weights
Wl		= rand(Dout, Nh).*w0*2-w0; %Output weights
Wnl2	= rand(Nh, Dout).*w0*2-w0; %Hidden weights
Wo		= rand(Ni, Nh).*w0*2-w0; %Output layer

J       = M;
gradJ   = 1e3;

while (gradJ > Theta),
    deltaWnl1	= 0;
    deltaWl		= 0;
    deltaWnl2	= 0;
    deltaWo		= 0;
    
    for m = 1:M,
        Xm = patterns(:,m);
        tk = patterns(:,m);
        
        %Forward propagate the input:
        %...to the first NL hidden units
        gh1				= Wnl1*Xm;
        [y1, dfh1]	= activation(gh1);
        
        %...to the linear hidden units
        y2				= Wl*y1;
        
        %...to the second NL hidden units
        gh3				= Wnl2*y2;
        [y3, dfh3]	= activation(gh3);

        %...to the output unit
        zk				= Wo*y3;;
        
        %Backpropagation!
        %Evaluate delta_k at the output: delta_k = (tk-zk)*f'(net)
        delta_k		= (tk - zk);
        
        %...and delta_j: delta_j = f'(net)*w_j*delta_k
        delta_j3		= dfh3.*sum((Wo.*(delta_k*ones(1,Nh))))';
        delta_j2		= sum((Wnl2.*(delta_j3*ones(1,Dout))))';
        delta_j1		= dfh1.*sum((Wl.*(delta_j2*ones(1,Nh))))';
        
        %delta_w_kj <- w_kj + eta*delta_k*y_j
        deltaWo		= deltaWo   + eta*delta_k*y3';
        deltaWnl2		= deltaWnl2 + eta*delta_j3*y2';
        deltaWl		= deltaWl   + eta*delta_j2*y1';
        deltaWnl1		= deltaWnl1 + eta*delta_j1*Xm';
        
    end
    
    %w_kj <- w_kj + eta*delta_Wo
    Wo				= Wo   + deltaWo;
    Wnl2				= Wnl2 + deltaWnl2;
    Wl				= Wl   + deltaWl;
    Wnl1				= Wnl1 + deltaWnl1;
    
    %Calculate total error
    oldJ = J;
    J    = 0;
    for i = 1:M,
       Xm				= patterns(:,i);
       gh1				= Wnl1*Xm;
       [y1, dfh1]		= activation(gh1);
       y2				= Wl*y1;
       gh3				= Wnl2*y2;
       [y3, dfh3]		= activation(gh3);
       zk				= Wo*y3;;
       J = J + sum((zk-Xm).^2);
    end
    J		 = sqrt(J/M);
    gradJ = abs(J - oldJ);
    
    iter 			= iter + 1;
    if (iter/IterDisp == floor(iter/IterDisp)),
        disp(['Iteration ' num2str(iter) ': Total error is ' num2str(J)])
    end
    
end

disp(['Backpropagation converged after ' num2str(iter) ' iterations.'])

%Translate the patterns by propagating them through the first two layers of the network
new_patterns = zeros(Dout,M);
for i = 1:M,
   Xm				= patterns(:,i);
   gh1			    = Wnl1*Xm;
   [y1, dfh1]	    = activation(gh1);
   new_patterns(:,i)= Wl*y1;
end




function [f, df] = activation(x)

a = 1.716;
b = 2/3;
f	= a*tanh(b*x);
df	= a*b*sech(b*x).^2;