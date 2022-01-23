%-------------------
% danny_bug1.m
% isolates bug in Danny's implementation.
% Dror, 3.31.2009
%-------------------

clear
%%% turn off randomness for debug
%rand('state',sum(100*clock));
%randn('state',sum(100*clock));
rand('state',0);
randn('state',0);

%-------------------
% Parameter definitions
%-------------------
num_reps=101;
iter=15; % Number of iterations in BP
n=500; %Signal length
k=n/10; %Sparsity
l=20; % # nonzeros per row
SNR=100;
sigma_1=sqrt(SNR); 
sigma_0=1;  %Noise IN the signal x
sigma_Z=1e-20;
epsilon=1e-10;
boundx=250;
model_order=243	% Number of sampling points 
R=15; % nonzeros per column

%-------------------
% Generate signal
%-------------------
disp ('GENERATING THE SIGNAL X...');
[x, heavyind]=generatex_noisy(n, k, sigma_1,sigma_0);
disp(sprintf('l2 norm of x: %g', norm(x) )); 

%-------------------
% Generate measurement matrix
%-------------------
[phi]=gen_phi(n, l,R);
phisign=randn(size(phi));  phisign=sign(phisign);
num_rows=size(phi,1);
phi2=zeros(num_rows,n);

 %-------------------
% Encode (compute measurements)
%-------------------
disp ('GENERATING THE MEASUREMENTS...');
phi2=zeros(size(phi,1),n); % direct representation of matrix
for ind_row=1:size(phi,1)
   phi2(ind_row,phi(ind_row,:))=phisign(ind_row,:);
end;
measvec=phi2*x'+sigma_Z*randn(size(phi2,1),1);	% add noise

%---------------
% compute priors and stuff
%---------------
[tmptmp, dispind]=sort(-abs(x)); % display stuff
dispind=[dispind(1:5), 1:5];

xx=1:model_order;   % ticks for computing pdf
xx=xx-(model_order+1)/2;  
xx=xx/max(xx)*boundx; % values over which pdf is sampled

pdf_prior= (k/n)*normpdf(xx,0,sigma_1); % signal prior
if (sigma_0 > epsilon)
  pdf_prior= pdf_prior  +    (1-k/n)*normpdf(xx,0,sigma_0);
else
  in2=find(abs(xx)<epsilon);
  pdf_prior(in2)=pdf_prior(in2) + (1-k/n);
end
pdf_prior=pdf_prior/sum(pdf_prior);

y_noise=normpdf(xx,0,sigma_Z);	% noise prior
y_noise=y_noise/sum(y_noise);

%---------------
% invoke Danny's code
%---------------
[xrecon, mrecon, srecon]=NBP(phi2,x,measvec,sigma_Z,iter,dispind,1e-10,pdf_prior,y_noise,xx);
