function Perror = Chernoff(mu1, sigma1, mu2, sigma2, p1)

% Find the Chernoff bound given means and covariances of single gaussian distributions
% Inputs:
%   mu1         - Mean for class 1
%   sigma1      - Covariance matrix for class 1
%   mu2         - Mean for class 2
%   sigma2      - Covariance matrix for class 2
%   p1          - Probability of class 1
%
% Outputs
%	Perror  	- Error bound

beta    = linspace(0,1,100);
k       = zeros(1,length(beta));

%First, claculate k(beta)
for i = 1:length(beta),
    k(i) = beta(i)*(1-beta(i))/2*(mu2-mu1)*inv(beta(i)*sigma1+(1-beta(i))*sigma2)*(mu2-mu1)'+...
              1/2*log(det(beta(i)*sigma1+(1-beta(i))*sigma2)/(det(sigma1)^beta(i)*det(sigma2)^(1-beta(i))));
end
      
%Find the minimum of exp(-k)
[m, index]  = min(exp(-k));
min_beta    = beta(index);

Perror      = p1^min_beta*(1-p1)^(1-min_beta)*exp(-k(index));