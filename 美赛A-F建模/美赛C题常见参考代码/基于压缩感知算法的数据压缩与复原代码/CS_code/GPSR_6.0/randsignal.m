
clc
clear
% N is the original signal length
N = 2^10;

% k is number of observations to make
M = 2^8;

% number of spikes to put down
% n_spikes = floor(.01*n);
n_spikes = 100;


% random +/- 1 signal
f = zeros(N,1);
q = randperm(N);
f(q(1:n_spikes)) = randn(n_spikes,1);
disp('Creating measurement matrix...');
A = randn(M,N);                         %π€≤‚æÿ’Û
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
plot(x)
