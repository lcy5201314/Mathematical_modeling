%�õ�һ֡������������ϡ���

function [hat_x,error]=voicedomp(Y,N)

K=1;
 Psi=fft(eye(N,N))/sqrt(N);    
 Ys=abs(Psi*Y);
 
 
% Ys=abs(fft(Y));
flag=max(Ys);
for ii=1:N
    
if  Ys(ii)>0.1
    K=K+1;
end
end
K
M=K*round(log2(N/K))+30;
m=2*K;
% �ָ��źź�ԭʼ�źŶԱ�
[hat_x,error]=omp(Y,M,N,m);


                      