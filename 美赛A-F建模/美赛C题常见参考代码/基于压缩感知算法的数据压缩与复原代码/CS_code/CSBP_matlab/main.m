%-------------------
% main.m
% Main file for Belief Propagation
%-------------------
% Original code by Shriram Sarvotham, 2006. 
% Cleaned up a bit by Dror Baron, November 2008.
%-------------------

clear
rand('state',sum(100*clock));

%-------------------
% Signal parameters
%-------------------
%Signal and noise parameters
n=200;  %Signal length
k=20;	%Sparsity
SNR=100;	% input snr
sigma_1=sqrt(SNR);
sigma_0=1; % small signal coefficients
sigma_Z=1; % noise in the measurements y (noisy measurements)

%-------------------
% CS-LDPC matrix
%-------------------
l=20; % constant row weight
r=10; % constant column weight

%-------------------
% CS-BP parameters
%-------------------
gamma_mdbpf=0.35; %Damping for Belief Prop
gamma_mdbpb=0.35;
gamma_pdbp=0.0;
iter=10; % Number of iterations in Belief Prop
p=243; % Number of sampling points (FFT runs decently for this value)

%-------------------
% Generate signal
%-------------------
t1=cputime;
disp ('GENERATING THE SIGNAL X...');
[x, heavyind]=generatex_noisy(n, k, sigma_1, sigma_0);
x=(x/norm(x))*sqrt(k);
x=sigma_1*x;
disp(sprintf('l2 norm of x: %g', norm(x) )); 
  
%-------------------
% Generate measurement matrix
%-------------------
[phi]=gen_phi(n, l, r);
phisign=randn(size(phi));  phisign=sign(phisign);

%-------------------
% Run driver, which decodes the signal
%-------------------
mrecon=driver_function(n,k,l,phi,phisign,x,SNR,sigma_1,sigma_0,sigma_Z,...
   iter,p,gamma_mdbpf,gamma_mdbpb,gamma_pdbp);
fprintf('error=%6.2f, time=%8.2f\n',norm(mrecon-x),cputime-t1);

plot(1:n,x,'-b',1:n,mrecon,'o');
xlabel('coefficient index');
ylabel('coefficient value');
legend('x','mrecon');