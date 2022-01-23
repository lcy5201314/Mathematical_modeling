function [mu, bias, varjack] = mean_bootstrap(data, B)

%Find the estimate of the mean, it's bias and variance using the bootstrap estimator method
%Inputs:
%	data	- The data from which to estimate
%	B		- Number of sets to draw
%Outputs:
%	mu		- The mean
%	bias	- The bias of the estimator
%	var	    - The variance of the estimate

[D, N] = size(data);

mu_star		= zeros(D,B);

for i = 1:B,
   %Draw N samples from the data, with replacement
   indices = zeros(1,N);
   for j = 1:N,
      indices(j) = 1 + floor(rand(1)*N);
   end
   
   %Average the data points chosen
   mu_star(:,i) = mean(data(:,indices)')';
end

mu			= 1/B*sum(mu_star')';						
bias		= mu-mean(data')';
varjack	= 1/B*sum((mu_star-mu*ones(1,B))'.^2)';