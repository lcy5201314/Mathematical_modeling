function [S1, S2] = Sufficient_statistics(data, distribution)

%Find sufficient statistics for data assuming a distribution
%Inputs:
%	data	- The data from which to estimate
%	distribution	- What distribution is the data taken from. Can be one of:
%							Normal, MultivariateNormal, Exponential, Rayleigh, Maxwell, Gamma, 
%							Beta, Poisson, Bernoulli, Binomial, Multinomial
%	
%Outputs:
%	S		- The sufficient statistics

[D, N] = size(data);

switch distribution
case 'Normal', 
   S1 = 1/N*sum(data')';
   S2 = 1/N*sum(data'.^2)';
case 'MultivariateNormal', 
   S1 = 1/N*sum(data')';
   S2 = 1/N*(data*data');
case {'Exponential', 'Poisson', 'Bernoulli', 'Binomial', 'Multinomial'}
   S1 = 1/N*sum(data')';
   S2 = [];
case {'Rayleigh', 'Maxwell'},
   S1	= 1/N*sum(data'.^2)';
   S2 = [];
case 'Gamma',
   S1	= prod(data').^(1/N);
   S2 = 1/N*sum(data')';
case 'Beta', 
   S1	= prod(data')'.^(1/N);
   S2 = prod(1-data')'.^(1/N);
otherwise
   error('Unknown distribution')
end
   