% romp vs omp

clc
clear
% N is the original signal length
N = 500;
% k is number of observations to make
M = 400;
% number of spikes to put down
% n_spikes = floor(.01*n);
k= 50;
% random +/- 1 signal
f = zeros(N,1);
q = randperm(N);
f(q(1:k)) = randn(k,1);
disp('Creating measurement matrix...');
A = randn(M,N);                         %π€≤‚æÿ’Û
T=orth(A')';
x=f;

y=T*x;



x1=zeros(N,1);
x2=zeros(N,1);
x3=zeros(N,1);

[r1,x1,k]=omp(T,N,y,k);
[r2,x2,iter2]=romp(k,T,y);
[r3,x3,iter3]=romp_orignal(k,T,y);



% disp('ploting...')
% subplot(311);plot(x1,'k.-');hold on;plot(x);title('omp')
% subplot(312);plot(x2,'k.-');hold on;plot(x);title('romp')
% subplot(313);plot(x3,'k.-');hold on;plot(x);title('romp_orignal')



error1=norm(x1-x)^2/norm(x)^2  ;
snr1=10*log10(1/error1)
error2=norm(x2-x)^2/norm(x)^2  ;
snr2=10*log10(1/error2)
error3=norm(x3-x)^2/norm(x)^2  ;
snr3=10*log10(1/error3)