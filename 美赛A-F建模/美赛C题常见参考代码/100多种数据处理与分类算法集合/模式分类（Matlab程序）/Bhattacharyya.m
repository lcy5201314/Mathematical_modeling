function Perror = Bhattacharyya(mu1, sigma1, mu2, sigma2, p1)

% Find the Bhattacharyya bound given means and covariances of single gaussian distributions
% Inputs:
%   mu1         - Mean for class 1
%   sigma1      - Covariance matrix for class 1
%   mu2         - Mean for class 2
%   sigma2      - Covariance matrix for class 2
%   p1          - Probability of class 1
%
% Outputs
%	Perror  	- Error bound

k_half = 1/8*(mu2-mu1)*inv((sigma1+sigma2)/2)*(mu2-mu1)'+1/2*log(det((sigma1+sigma2)/2)/sqrt(det(sigma1)*det(sigma2)));

Perror = sqrt(p1*(1-p1))*exp(-k_half);

