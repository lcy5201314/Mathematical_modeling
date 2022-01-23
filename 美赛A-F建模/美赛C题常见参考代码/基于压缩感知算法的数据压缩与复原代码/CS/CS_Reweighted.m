%NAME: Reweighted L1-Minimization Tester
%PURPOSE: Tests sparse and noisy signals on Reweighted L1
%AUTHOR: Deanna Needell
%OUTSIDE FUNCTIONS: CVX package (Michael Grant and Stephen Boyd)
N = 256; %dimension
M = 128; %measurements
kVals = [30]; %sparsity
eepsVals = [1];
numTrials = 1;
maxIter = 9;
errorVecDecoding = zeros(length(kVals),length(eepsVals),maxIter,numTrials);
errors = zeros(numTrials,1);
for trial = 1:numTrials
for kIndex = 1:length(kVals)
K = kVals(kIndex);
for eIndex = 1:length(eepsVals)
eeps = eepsVals(eIndex);
% Gaussian spikes in random locations
x = zeros(N,1); q = randperm(N);
x(q(1:K)) = sign(randn(K,1));
% measurement matrix
Phi = sign(randn(M,N))/sqrt(M);
% observations
err = randn(M, 1);
sigma = 0.2*norm(Phi*x,2)/norm(err,2);
err = sigma*err;
y = Phi*x + err;

errors(trial) = norm(err,2);




for iter = 1:maxIter

    
%Set the weights
if iter > 1
weights = 1./(abs(xDecoding)+eeps/(1000*iter));
else
weights = 1*ones(N,1);
end


%Set noise tolerance parameter
delta=sqrt(sigma^2*(M+2*sqrt(2*M)));
%Use CVX to perform minimization
cvx_begin
cvx_quiet(true)
variable xa(N);
minimize( norm(diag(weights)*xa,2) );
subject to
norm(Phi*xa - y, 2) <= delta;
cvx_end
%Collect results
xDecoding = xa;
errorVecDecoding(kIndex,eIndex,iter,trial) = norm((x-xDecoding),2);
end


plot(x,'k-');
hold on;
plot(xDecoding,'r-');


end
end
end

