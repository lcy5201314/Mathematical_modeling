function a = stochastic_regression(train_patterns, train_targets, param)

% Find the linear regression using stochastic regression
% Inputs:
% 	patterns	- Train patterns
%	targets	    - Train targets
%	param		- [Maximum iteration, Theta (Convergence criterion)]
%
% Outputs
%	a			- Line equation

[c, n]		   = size(train_patterns);

Max_iter	   = param(1);
theta		   = param(2);

%Initial weights
a              = zeros(c,1);
R              = inv(train_patterns*train_patterns');
iter  	       = 0;

update	       = 1e3;

y              = train_patterns;

while ((sum(abs(update)) > theta) & (iter < Max_iter))
   iter = iter + 1;
   
   %Randomally choose a sample
   i    = randperm(n);
   k    = i(1);

   %R_k+1 = R_k + R_k*(R_k*y_k)'/(1+y_k'*R_k*y_k)
   R    = R - R*y(:,k)*(R*y(:,k))';
   
   % a <- a + R_k+1*(Target(k) - a'*y_k)*y_k
   delta= R*(train_targets(k) - a'*y(:,k))*y(:,k);
   a    = a + delta;
   
   update = sum(abs(delta));
end

if (iter == Max_iter),
   disp(['Maximum iteration (' num2str(Max_iter) ') reached']);
else
   disp(['Did ' num2str(iter) ' iterations'])
end

