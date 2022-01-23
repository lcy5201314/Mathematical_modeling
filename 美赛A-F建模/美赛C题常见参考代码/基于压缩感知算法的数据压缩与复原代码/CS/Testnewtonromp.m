tic
clc
clear
% N is the original signal length
N = 2^10;
% k is number of observations to make
M = 450;

% number of spikes to put down
% n_spikes = floor(.01*n);
k= 120;


% random +/- 1 signal
f = zeros(N,1);
q = randperm(N);
f(q(1:k)) = randn(k,1);
disp('Creating measurement matrix...');
A = randn(M,N);                         %π€≤‚æÿ’Û
T=orth(A')';

% for kk=2:N
%     for nn=1:N
%         dctbasis(kk,nn)=(2/N)^0.5*cos((2*(nn-1)+1)*(kk-1)*pi/2/N);
%     end
% end
% for nn=1:N
%     dctbasis(1,nn)=(1/N)^0.5*cos((2*(nn-1)+1)*(1-1)*pi/2/N);
% end
% B=dctbasis;                        %  ∏µ¿Ô“∂’˝±‰ªªæÿ’Û
x=f;

y=T*x;
figure(1)
plot(x,'r');
hold on;

Aug_t=[];
rn=y;

hat_x=zeros(N,1);
lamada=1/2;

[r,hat_x,numIts] = newtonromp(k, T, y,lamada);

norm(rn)
numIts
error=norm(hat_x-x)^2/norm(x)^2 ;
snr=10*log10(1/error)

disp('draw the recovery signal ...')
plot(hat_x,'k.-');
hold off
toc