% clc;
% clear;
% close
% %  1. 时域测试信号生成
% K=16;       %  稀疏度(做FFT可以看出来)
% N=320;      %  信号长度
% M=100;      %  测量数(M>=K*log(N/K),至少40,但有出错的概率)
% m=2;%迭代次数
% % f1=50;    %  信号频率1
% % f2=100;   %  信号频率2
% % f3=200;   %  信号频率3
% % f4=400;   %  信号频率4
% % fs=800;   %  采样频率
% % ts=1/fs;  %  采样间隔
% % Ts=1:N;   %  采样序列
% % x=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.1*cos(2*pi*f3*Ts*ts)+0.9*cos(2*pi*f4*Ts*ts);  %  完整信号
% % x=x';
% [x,fs,bits]=wavread('M1_1',[8400 8719]);
% %  2.  时域信号压缩传感
% for kk=2:N
%     for nn=1:N
%         dctbasis(kk,nn)=(2/N)^0.5*cos((2*(nn-1)+1)*(kk-1)*pi/2/N);
%     end
% end
% for nn=1:N
%     dctbasis(1,nn)=(1/N)^0.5*cos((2*(nn-1)+1)*(1-1)*pi/2/N);
% end
% 
% Psi=dctbasis;
%    
% Phi=randn(M,N);  
% T=Phi*Psi'; 
tic%  测量矩阵(高斯分布白噪声)
clc
clear
% N is the original signal length
N = 2^10;
% k is number of observations to make
M = 450;

% number of spikes to put down
% n_spikes = floor(.01*n);
k= 150;


% random +/- 1 signal
f = zeros(N,1);
q = randperm(N);
f(q(1:k)) = randn(k,1);
disp('Creating measurement matrix...');
A = randn(M,N);                         %观测矩阵
T=orth(A')';
x=f;
y=T*x;                                          %  获得线性测量 
m=5;
[xr iter_num res]=SP(y,T,k,m);
% x_recovery=real(Psi'*xsp);
rn=norm(res)
iter_num
plot(x,'r');

hold on;
plot(xr','k.-')
legend('Original','Recovery')
error=norm(xr-x)^2/norm(x)^2  ;
snr=10*log10(1/error)
hold off
toc

