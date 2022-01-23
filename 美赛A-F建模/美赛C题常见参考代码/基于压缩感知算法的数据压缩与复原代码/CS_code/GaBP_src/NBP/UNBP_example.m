%Code written by Danny Bickson
%Code is found on: http://www.cs.huji.ac.il/labs/danss/p2p/gabp/index.html
%Example of solving boolean least squares via NBP
%The problem is argmin_x ( norm (Hx - y))  s.t. x \in {-1,1}
% G is 121 x 121 matrix. Ulike ldlc G is the sparse ENCODING matrix (and not the decoding)

function [] = UNBP_example()
clear
rand('state',sum(100*clock));
randn('state',sum(100*clock));


load('H_n_121_d_3_1.mat','H');
G = H;
n=size(H,1);

%addpath('../CSBP_matlab');
max_rounds = 12; % max number of iterations
epsilon = 1e-70; % 
transmit = rand(n,1);
transmit(transmit<0.5) = -1;
transmit(transmit>=0.5) = 1;
sigma = 0.1;

y = G*transmit+randn(n,1)*sigma;


epsilon=1e-70;
boundx=5;
%---------------
% compute signal node prior
%---------------
model_order = 1543;
xx=1:model_order;  xx=xx-(model_order+1)/2;  xx=xx/max(xx);
xx=xx*boundx; % values over which pdf is sampled

pdf_prior= normpdf(xx,-1,.2) + normpdf(xx,1,.2);
%pdf_prior= normpdf(xx,0,1);
pdf_prior=pdf_prior/sum(pdf_prior);


disp ('STARTING THE DECODER...');

displayind = [ 1 2 3];

self_pot=cell(1,2*n);
BIGG=[sparse(n,n) G; G' sparse(n,n)];
for i=1:n
   self_pot{i} = pdf_prior; 
   self_pot{i+n} = normpdf(xx,y(i),sigma);
end


 [xrecon, mrecon, srecon]=UNBP(BIGG,transmit,y,self_pot,...
    max_rounds, displayind ,epsilon, xx);


fprintf('[   %d   %1d   %d   %d]\n', transmit(displayind));
fprintf('success=%6.2f \n',sum((((mrecon(1:n)>0)*2)-1)==transmit')/n);

fprintf('GaBP success=%6.2f  \n',sum(transmit == sign(inv(G)*y))/n);


end