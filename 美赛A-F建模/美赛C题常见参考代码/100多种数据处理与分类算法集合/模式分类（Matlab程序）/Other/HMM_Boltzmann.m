function [a, b] = HMM_Boltzmann(Nh, No, V, init_state)

% Find the probability transition matrices a,b from sample data using the Boltzmann network algorithm
%
% Inputs:
%	Nh					- Number of hidden states
%	No					- Number of output states 
%	V					- Observed output sequence
%  init_state     - Initial state of the HMM
%
% Output:
%	a					- Estimated transition probability matrix 
%	b					- Estimated output generator matrix 
%
%	In this implementation we hold two separate state matrices: Sh - The hidden state matrix
%	and So - The output state matrix. This is possible because there are no direct 
%	connections between the input (previous time step) and output (Visible state) of the chain

theta       = 1e-6;
eta_anneal  = 0.95;
Tini        = 5;
Tmin        = 0.1;
Tf				= length(V);
Tb				= Tini;
iter			= 0;
max_iter		= 1e5;
DispIter		= 10;

%Make symmetric weight matrices with a zero diagonal
a				= rand(Nh);
b				= rand(Nh, No);
W				= [zeros(Nh), a, zeros(Nh, No); a', zeros(Nh), b; zeros(No, Nh), b', zeros(No,No)];
zero_diag	= ~eye(size(W));

while ((iter < max_iter) & (Tb > Tmin)),
    
    %Randomize states Si
    %The states are the hidden state values. Therefore, only one value per time step can
    %be active
    Sh	= rand(Nh, Tf);
    So	= rand(No, Tf);
    for i = 1:Tf,
       Sh(:,i) = (Sh(:,i) == max(Sh(:,i)));
       So(:,i) = (So(:,i) == max(So(:,i)));
    end
    
    %Clamp outputs
    for i = 1:Tf,
       So(:,i) 	= zeros(No, 1);
       So(V(i),i) = 1;
    end
    
    Sin	= [Sh; So];
    
    %Anneal network with input and output clamped
    T = Tini;
    while (T > Tmin),
       for i = 1:Tf,
          %Select a node randomally
          j	= floor(rand(1)*Nh)+1;
          
          %li<-sigma(w_ij*s_j)
          if (i == 1),
             x					= zeros(Nh, 1);
             x(init_state)	= 1;
          else
             x					= Sin(1:Nh,i-1);
          end
	       S		= [x; Sin(:,i)];   
          l		= W(:,j+Nh)'*S;
          
          %Si<-f(li, T(k))
          Sin(j,i) = tanh(l/T);
       end
       
       %Lower the temperature
       T = T*eta_anneal;   
    end
    
    %At final T calculate [SiSj]alpha_i,alpha_o clamped
    for i = 1:Tf,
       if (i == 1),
          x					= zeros(Nh, 1);
          x(init_state)	= 1;
       else
          x					= Sin(1:Nh,i-1);
       end
       S								= [x; Sin(:,i)];   
       SiSj_io_clamped(i,:,:) = S*S';
    end
    
    %Randomize states Si
    Sh	= rand(Nh, Tf);
    So	= rand(No, Tf);
    for i = 1:Tf,
       Sh(:,i) = (Sh(:,i) == max(Sh(:,i)));
       So(:,i) = (So(:,i) == max(So(:,i)));
    end
    Sin	= [Sh; So];
    
    %Anneal network with input clamped and output free
    T = Tini;
    while (T > Tmin),
       for i = 1:Tf,
          %Select a node randomally
          j	= floor(rand(1)*(Nh+No))+1;
          
          %li<-sigma(w_ij*s_j)
          if (i == 1),
             x					= zeros(Nh, 1);
             x(init_state)	= 1;
          else
             x					= Sin(1:Nh,i-1);
          end
	       S		= [x; Sin(:,i)];   
          l		= W(:,j+Nh)'*S;
          
          %Si<-f(li, T(k))
          Sin(j,i) = tanh(l/T);
       end
       
       %Lower the temperature
       T = T*eta_anneal;   
    end
    
    %At final T calculate [SiSj]alpha_i clamped
    for i = 1:Tf,
       if (i == 1),
          x					= zeros(Nh, 1);
          x(init_state)	= 1;
       else
          x					= Sin(1:Nh,i-1);
       end
       S								= [x; Sin(:,i)];   
       SiSj_i_clamped(i,:,:)  = S*S';
    end
    
    %Update W
    %Wij<-Wij + eta/T([SiSj]alpha_i,alpha_o clamped - [SiSj]alpha_i clamped)
    dW = squeeze(mean(SiSj_io_clamped - SiSj_i_clamped)).*zero_diag;
    W  = W + eta_anneal*Tb*dW;
    
    Tb 	= Tb * eta_anneal;
    iter	= iter + 1;
    if (floor(iter/DispIter) == iter/DispIter),
       disp(['Iteration ' num2str(iter) ': Temperature is ' num2str(Tb)])
    end
    
end
 
A 	= W(1:Nh, Nh+1:2*Nh);
B	= W(Nh+1:2*Nh, 2*Nh:end);
 
a	= exp(A/T);
b	= exp(B/T);

a	= a ./ (sum(a')' * ones(1,size(a,2)));
b	= b ./ (sum(b')' * ones(1,size(b,2)));
