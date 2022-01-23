%Code written by Danny Bickson
%Code is found on: http://www.cs.huji.ac.il/labs/danss/p2p/gabp/index.html
%Example of solving boolean least squares via NBP
%The problem is argmin_x ( norm (Gx - y))  s.t. x \in {-1,1}
%G is a 3x3 matrix
% Note that unlike LDLC paper by Sommer et al., G is a sparse ENCODING matrix

function [] = UNBP_small_example()
clear
rand('state',sum(100*clock));
randn('state',sum(100*clock));


% G=[0 0.5 0.8; 0.5 0 -0.5; 0.8 -0.5 0];
% n=size(G,1);
% 
% %addpath('../CSBP_matlab');
% max_rounds = 12; % max number of iterations
% epsilon = 1e-70; % 
% transmit = rand(n,1);
% transmit(transmit<0.5) = -1;
% transmit(transmit>=0.5) = 1;
% sigma = 1;
% 
% y = G*transmit+randn(n,1)*sigma;


epsilon=1e-70;
boundx=5;
%---------------
% compute signal node prior
%---------------
model_order = 1111;
xx=1:model_order;  xx=xx-(model_order+1)/2;  xx=xx/max(xx);
xx=xx*boundx; % values over which pdf is sampled

pdf_prior= normpdf(xx,-1,.2) + normpdf(xx,1,.2);
%pdf_prior= normpdf(xx,0,1);
pdf_prior=pdf_prior/sum(pdf_prior);




displayind = [ 1 2 3];


max_rounds = 20; % max number of iterations
%epsilon = 0.001; % convergence threshold

H=[0 0.4 0.5; 0.4 0 -0.4; 0.5 -0.4 0];
% BIGH = [ eye(3,3) H; H' eye(3,3)];
y =ones(3,1);
% BIGy=[zeros(3,1)' y']';

 [x,p] = GBP(H+eye(3),y,3,0.00000001);
% x
% p
%self_pot_val = kde(1,1);
self_pot_val=normpdf(xx,0,1);
self_pot_val=self_pot_val/sum(self_pot_val);
self_pot=cell(1,3);
self_pot{1} = self_pot_val;
self_pot{2} = self_pot_val;
self_pot{3} = self_pot_val;
disp ('STARTING THE DECODER...');
%x = NBP(H,self_pot,30, 1e-10,'exact')
 [meanx, varx, maxx,X]=UNBP(H,(inv(H)*y),y,self_pot,...
    3, displayind ,epsilon, xx);

%comparing the result to iterative linear detection. See:
%Gaussian belief propagation based multiuser detection. D. Bickson, O.
%Shental, P. H. Siegel, J. K. Wolf, and D. Dolev, In IEEE Int. Symp. on Inform. 
%Theory (ISIT), Toronto, Canada, July 2008. 

disp(['answer should be']);
x
disp(['anser is ']);
meanx
varx
maxx

% self_pot=cell(1,2*n);
% BIGG=[sparse(n,n) G; G' sparse(n,n)];
% for i=1:n
%    self_pot{i} = pdf_prior; 
%    self_pot{i+n} = normpdf(xx,y(i),sigma);
% end
% 
% 
%  [xrecon, mrecon, srecon,X]=UNBP(BIGG,transmit,y,self_pot,...
%     max_rounds, displayind ,epsilon, xx);


% fprintf('[   %d   %1d   %d   %d]\n', transmit(displayind));
% fprintf('success=%6.2f \n',sum((((mrecon(1:n)>0)*2)-1)==transmit')/n);
% 
% fprintf('GaBP success=%6.2f  \n',sum(transmit == sign(inv(G)*y))/n);


end