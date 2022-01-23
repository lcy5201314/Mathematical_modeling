function out = ALM_SADAL_smoothed(D, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a sample code for testing the algorithm in our paper 
% ''Fast Alternating Linearization Methods for
%       Minimizing the Sum of Two Convex Functions'', Donald Goldfarb,
%       Shiqian Ma and Katya Scheinberg, Tech. Report, Columbia University,
%       2009 - 2010. Preprint available at: 
%       http://arxiv.org/pdf/0912.4571v2.pdf
%
% Author: Shiqian Ma
% Date  : Apr. 20, 2010 
% IEOR, Columbia University, Copyright (2010)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n] = size(D);
mu = opts.mu; sigma = opts.sigma; rho = opts.rho;
Y = zeros(m,n); gradgY = zeros(m,n);
X = D - Y; Dnorm = norm(D,'fro');
Lambda = zeros(m,n); sv = opts.sv;

for itr = 1: opts.maxitr
    [U, gamma, V] = svd(mu*Lambda-Y+D, 'econ');
    gamma = diag(gamma);
    gamma_new = gamma-mu*gamma./max(gamma,mu+sigma);
    
    X = U*diag(gamma_new)*V';
    Lambda = Lambda - (X+Y-D)/mu;
    muY = mu;
    B = Lambda - (X-D)/muY;
    
    Y = muY*B - muY*min(rho, max(-rho,muY*B/(sigma+muY)));
    Lambda = Lambda - (X+Y-D)/muY;

    StopCrit = norm(D-X-Y,'fro')/Dnorm;
    fprintf('iter: %d, crit:%3.2e\n', ...
        itr, norm(D-X-Y,'fro')/norm(D,'fro'));

    mu = max(opts.muf, mu*opts.eta_mu);
    sigma = max(opts.sigmaf, sigma*opts.eta_sigma);


    if StopCrit < opts.epsilon
        out.X = X; out.Y = Y; out.iter = itr; out.StopCrit = StopCrit;
        return;
    end
end
out.X = X; out.Y = Y; out.iter = itr; out.StopCrit = StopCrit;
