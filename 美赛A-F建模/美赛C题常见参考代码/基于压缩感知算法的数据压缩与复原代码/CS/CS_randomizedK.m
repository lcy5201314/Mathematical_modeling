%NAME: Randomized Kaczmarz Tester
%PURPOSE: Tests RK on noisy systems
%AUTHOR: Deanna Needell
%OUTSIDE FUNCTIONS: none
clear all
warning off all
m = 100; %rows
n=100; %columns
numIts = 1000;
numTrials = 100;
A = zeros(numTrials, m, n);
e = zeros(numTrials, m);
x = zeros(numTrials, n);
b = zeros(numTrials, m);
est = zeros(numTrials, numIts, n); %estimations
initErr = zeros(numTrials);
R = zeros(numTrials, 1); %Value of R (as in paper)
mu = zeros(numTrials); %Coherence
gamma = zeros(numTrials, 1); %worst error to row norm ratio
beta = zeros(numTrials); %beta is in theorem
errorsRK = zeros(numTrials, numIts);
errorsCG = zeros(numTrials, numIts);
for trial=1:numTrials
if mod(trial, 1) == 0
display(['Trial ', num2str(trial)]);
end
%Set matrix equation
A(trial, :, :) = sign(randn(m,n));
x = zeros(n, 1)';
b(trial,:) = reshape(A(trial, :, :), m, n)*(x');
%Set initial guess, ||x - x0|| = 1
est(trial, 1, :) = randn(1, n);
est(trial, 1, :) = est(trial, 1, :) / norm(reshape(est(trial, 1, :), 1, n),2);
origest = reshape(est(trial, 1, :), 1, n)';
%Add error to RHS of Ax=0
e(trial, :) = randn(1, m)*2;
e(trial, :) = e(trial, :) / norm(e(trial, :), 2) / 10;
b(trial, :) = b(trial, :)+e(trial, :); %%%%NOISY!!
%Calculate stats
initErr = norm(reshape(est(trial, 1, :), 1, n) - x,2);
fronorm = norm(reshape(A(trial, :, :), m, n), 'fro');
R(trial) = norm(pinv(reshape(A(trial, :, :), m, n)),2)*fronorm;
temp = zeros(1, m);
for i=1:m
temp(i) = norm(reshape(A(trial, i, :), 1, n), 2);
end
gamma(trial) = max(abs( e(trial, :)./temp ));
errorsRK(trial, 1) = norm(reshape(est(trial, 1, :), 1, n)-x,2);
errorsCG(trial, 1) = norm(reshape(est(trial, 1, :), 1, n)-x,2);
%Run RKA.6. Randomized Kaczmarz 162
for it=2:numIts
%Select random hyperplane
pick = rand * fronorm^2;
counter = 0;
index = 1;
while counter + norm(reshape(A(trial, index, :),1, n), 2)^2 < pick
counter = counter + norm(reshape(A(trial, index, :),1, n), 2)^2;
index = index + 1;
end
%Modify estimation
est(trial, it, :) = est(trial, it-1, :) + (b(index) - dot((A(trial, index,est(trial, it-1, :)))/ norm(reshape(A(trial, index, :), 1, n),2)^2));
errorsRK(trial, it) = norm(reshape(est(trial, it, :), 1, n)-x,2);
end
%Run CG
for it=1:numIts
[estcg,flag] = cgs(reshape(A(trial, :, :),m, n),b(trial,:)',10^(-10), it)
errorsCG(trial, it) = norm(estcg-x',2);
end
end