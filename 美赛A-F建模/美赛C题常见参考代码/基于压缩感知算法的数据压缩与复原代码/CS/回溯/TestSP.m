clc;
clear;

%  1. 时域测试信号生成
K=30;       %  稀疏度(做FFT可以看出来)
N=320;      %  信号长度
M=100;      %  测量数(M>=K*log(N/K),至少40,但有出错的概率)
m=18;
% f1=50;    %  信号频率1
% f2=100;   %  信号频率2
% f3=200;   %  信号频率3
% f4=400;   %  信号频率4
% fs=800;   %  采样频率
% ts=1/fs;  %  采样间隔
% Ts=1:N;   %  采样序列
% x=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.1*cos(2*pi*f3*Ts*ts)+0.9*cos(2*pi*f4*Ts*ts);  %  完整信号
% x=x';
[x,fs,bits]=wavread('M1_1',[8400 8719]);
%  2.  时域信号压缩传感
Psi=fft(eye(N,N))/sqrt(N);                        %  傅里叶正变换矩阵
   
Phi=randn(M,N);  
T=Phi*Psi';                                       %  测量矩阵(高斯分布白噪声)
y=Phi*x;                                          %  获得线性测量 
sigma=1;
[xsp iter_num res]=SP(y,T,K,m);
x_recovery=real(Psi'*xsp);
rn=norm(res)
iter_num
figure(2)
plot(x_recovery','k.-')
hold on;
plot(x,'r');
legend('Recovery','Original')
error=norm(x_recovery-x)^2/norm(x)^2  ;
snr=10*log10(1/error)
hold off


