function [mu, bias, varjack] = mean_jackknife(data)

%Find the estimate of the mean, it's bias and variance using the jackknife estimator method
%Inputs:
%	data	- The data from which to estimate
%Outputs:
%	mu		- The mean
%	bias	- The bias of the estimator
%	var	    - The variance of the estimate

[D, N] = size(data);

mu_i		= zeros(D,N);

for i = 1:N,
   mu_i(:,i) = mean(data(:,[1:i-1,i+1:N])')';
end

mu			= 1/N*sum(mu_i')';
bias		= (N-1)*(mu-mean(data')');
varjack	= (N-1)/N*sum((mu_i-mu*ones(1,N))'.^2)';