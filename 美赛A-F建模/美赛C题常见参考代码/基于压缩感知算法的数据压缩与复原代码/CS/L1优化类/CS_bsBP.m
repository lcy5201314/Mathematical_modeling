%NAME: Basis Pursuit Tester²©Ê¿ÂÛÎÄ
%PURPOSE: Tests sparse, noisy, compressible signals on Basis Pursuit
%AUTHOR: Deanna Needell
%OUTSIDE FUNCTIONS: l1eq_pd (L1-Magic, J. Romberg)

clc;
clear all
warning off all
%Variables
N=[10:5:250];
d=[256];
n=[2:2:90];
p=[0.2:0.2:1]; %Used for compressible signals
%Number of Trials for each parameter set
numTrials = 500;

%Counters
numN = size(N, 2);
numd = size(d, 2);
numn = size(n, 2);
nump = size(p, 2);
%Data Collection
numCorr = zeros(numN, numd, numn);
minMeas = zeros(numN);
AvgerrorNormd = zeros(numN, numd, numn, nump);
Avgerror = zeros(numN, numd, numn, nump);
for ip=1:nump 
for id = 1:numd
for iN=1:numN
in=1;
done=0;
while in <= numn && ~done
for trial = 1:numTrials
tN = N(1, iN);
td = d(1, id);
tn = n(1, in);
tp = p(1, ip);
%Set Matrix
Phi = randn(tN, td);
Phi = sign(Phi);
I = zeros(1,1);
%Set signal
z = randperm(td);
v = zeros(td, 1);
for t = 1:tn
v(z(t))=1;
end
%USED IN THE CASE OF COMPRESSIBLE SIGNALS
%y = sign(randn(td, 1));
%noiErr=0;
%for i=1:tn
%v(z(i)) = i^(-1/tp)*y(i);
%if i > tn
% noiErr = noiErr + abs(v(z(i)));
%end
%end
%Set measurement and residual
x = Phi* v;
%USED IN THE CASE OF NOISE
%e = randn(tN, 1);
%x = x + e / norm(e, 2) / 2;
%Set initial estimate
x0 =Phi'*x;
%Run Basis Pursuit (via L1-Magic)
xp = l1eq_pd(x0, Phi, [], x, 1e-6);
%Collect Data
if norm(xp-v, 2) < 0.01
numCorr(iN, id, in) = numCorr(iN, id, in) +1;
end
AvgerrorNormd(iN, id, in, ip) = (AvgerrorNormd(iN, id, in, ip) *(counted-1) + (norm(xp-v,2)/noiErr*(tn)^0.5))/counted;
Avgerror(iN, id, in, ip) = (Avgerror(iN, id, in, ip) *(counted-1) + norm(xp-v, 2))/counted;
end % end Trial
if numCorr(iN, id, in) / numTrials > 0.98
minMeas(iN) = tn;
else
done=1;
end
in = in +1;
end % n
end % N
end % d
end % p