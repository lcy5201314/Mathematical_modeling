function [test_targets, updates] = Deterministic_Boltzmann(train_patterns, train_targets, test_patterns, params);

% Classify using the deterministic Boltzmann algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	Params  - [Ni, Nh, eta, Type, param], where:
%						Ni 	- Number of input units
%						Nh 	- Number of hidden units
%						eta	- Cooling rate
%						Type	- Type of weak learner
%                       params  - Weak learner parameters
%
% Outputs
%	test_targets	- Predicted targets
%   updates         - The updates to the network through the training iterations
%
% In this implementation, we train Ni weak learners, and find the best combination of 
% these classifiers

[Nio,M]                     = size(train_patterns);
No			                = 1;
[Ni, Nh, eta, type, Wparams]= process_params(params);

eta_anneal  = 0.95;	%Eta for the anneal process for the examples
iter		= 0;
max_iter    = 1e5;
DispIter    = 50;
update      = 0;
err		    = 1;
min_err	    = 0.01;

%First, build the weak learners and find the appropriate label for each input pattern 
%for these classifiers
weak_test_targets   = zeros(Ni, size(test_patterns,2));
weak_train_targets  = zeros(Ni, M);
for i = 1:Ni,
    weak_targets            = feval(type, train_patterns, train_targets, [train_patterns, test_patterns], Wparams);
    weak_train_targets(i,:) = weak_targets(1:M);
    weak_test_targets(i,:)  = weak_targets(M+1:end);
end

patterns = weak_train_targets*2-1;	%Inputs have to be {-1,+1}
targets  = train_targets*2-1;

%BEGIN BOLTZMANN LEARNING
Tini	    = max(eig(cov(patterns',1)'))/2;    %Initial temperature
Tmin        = 0.01;                             %Stopping temperature
Tb		    = Tini;

%Make a symmetric weight matrix with a zero diagonal
zero_diag   = ones(Ni+Nh+No) - eye(Ni+Nh+No);
W		    = rand(Ni+Nh+No)*2*sqrt(3/(Ni+Nh+No))-sqrt(3/(Ni+Nh+No));
W		    = W.*zero_diag;
for i = 1:Ni+Nh+No,
    for j = i:Ni+Nh+No,
        W(i,j) = W(j,i);
    end
end


while ((iter < max_iter) & (Tb > Tmin)),
   %Randomally select a training pattern
   in	= floor(rand(1)*M)+1;
   x	= patterns(:,in);
   t  = targets(in);
   
   %Randomize input states Si
   Sin = (rand(Nh,1)>0.5)*2-1;
   
   %Anneal network with input and output clamped
   T = Tini;
   while (T > Tmin),
      %Select a node randomally
      i	= floor(rand(1)*Nh)+1;
      
      %li<-sigma(w_ij*s_j)
      S	= [x; Sin; t]; 
      l	= W(:,i+Ni)'*S;
      
      %Si<-f(li, T(k))
      Sin(i)	= tanh(l/T);
      
      %Lower the temperature
      T = T*eta_anneal;   
   end
   
   %At final, low T, calculate [SiSj]alpha_i,alpha_o clamped
   S	= [x; Sin; t]; 
   SiSj_io_clamped = S*S';
   
   %Randomize input states Si
   Sin= (rand(Nh+No,1)>0.5)*2-1;
   
   %Anneal network with input clamped but output free
   T = Tini;
   while (T > Tmin),
      %Select a node randomally
      i	= floor(rand(1)*(Nh+No))+1;
      
      %li<-sigma(w_ij*s_j)
      S	= [x; Sin]; 
      l	= W(:,i+Ni)'*S;
      
      %Si<-f(li, T(k))
      Sin(i)	= tanh(l/T);
      
      %Lower the temperature
      T = T*eta_anneal;   
   end
   
   %At final, low T, calculate [SiSj]alpha_i clamped
   S	= [x; Sin];
   SiSj_i_clamped = S*S';
   
   %Wij<-Wij + eta/T([SiSj]alpha_i,alpha_o clamped - [SiSj]alpha_i clamped)
   dW = (SiSj_io_clamped - SiSj_i_clamped).*zero_diag;
   W  = W + eta*Tb*dW;
   
   Tb 	= Tb * eta;
   iter	= iter + 1;
   updates(iter) = sum(sum(abs(dW)))/size(dW,1)^2;
   if (floor(iter/DispIter) == iter/DispIter),
      disp(['Iteration ' num2str(iter) ': Temperature is ' num2str(Tb), ' Average weight update is: ' num2str(update/DispIter)])
      update = 0;
  else
      update = update + updates(iter);
   end
   
end

%Classify test patterns 
weak_test_targets = 2*weak_test_targets - 1;

%Randomize input states Si
Sin= (rand(Nh+No,size(test_patterns,2))>0.5)*2-1;

%Anneal network with input clamped but output free
T = Tini;
while (T > Tmin),
   %Select a node randomally
   i	= floor(rand(1)*(Nh+No))+1;
   
   %li<-sigma(w_ij*s_j)
   S	= [weak_test_targets; Sin]; 
   l	= W(:,i+Ni)'*S;
   
   %Si<-f(li, T(k))
   Sin(i,:)	= tanh(l/T);
   
   %Lower the temperature
   T = T*eta_anneal;   
end

So = Sin(Nh+1:end,:) > 0;

test_targets = So;