%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%重构概率newton,omp,romp,r_newton
clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%genarate signal
% N is the original signal length
N = 500;
% k is number of observations to make
M = 200;
% number of spikes to put down
% n_spikes = floor(.01*n);
k= 50;
% random +/- 1 signal
x = zeros(N,1);
q = randperm(N);
x(q(1:k)) = randn(k,1);
disp('Creating measurement matrix...');


x1=zeros(N,1);
x2=zeros(N,1);
x3=zeros(N,1);

n=20;
n_1=[];
n_2=[];
n_3=[];
rate=0.1:0.05:0.5;
for rr=0.1:0.05:0.5
    M=rr*N;
    A = randn(M,N);                         %观测矩阵
    T=orth(A')';
   
    y=T*x;
    n1=0;
    n2=0;
    n3=0;
    for jj=1:n
        [r1,x1,k]=omp(T,N,y,k);
        [r2,x2,iter2]=romp(k,T,y);
        [r3,x3,iter3]=romp_orignal(k,T,y);
        
        
        error1=norm(x1-x)^2/norm(x)^2  ;
        snr1=10*log10(1/error1);
        error2=norm(x2-x)^2/norm(x)^2  ;
        snr2=10*log10(1/error2);
        error3=norm(x3-x)^2/norm(x)^2  ;
        snr3=10*log10(1/error3);
        
        
        if snr1>100
            n1=n1+1;
        end
        if snr2>100
            n2=n2+1;
        end
        if snr3>100
            n3=n3+1;
        end
        
        
    end
    n_1=[n_1,n1/n];
    n_2=[n_2,n2/n];
    n_3=[n_3,n3/n];
    
end

plot(rate,n_1,'r*-');hold on;
plot(rate,n_2,'ko-')
plot(rate,n_3,'b.-')
ylabel('Recovery percentage')
xlabel('rate')
legend('omp','romp','romp-orignal')
