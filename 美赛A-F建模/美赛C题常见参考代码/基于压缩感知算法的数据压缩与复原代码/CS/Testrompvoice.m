clc;
clear;

%%  1. 时域测试信号生成
K=10;      %  稀疏度(做FFT可以看出来)
N=320;    %  信号长度
M=100;     %  测量数(M>=K*log(N/K),至少40,但有出错的概率)
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
%%  2.  时域信号压缩传感
Phi=randn(M,N);                                   %  测量矩阵(高斯分布白噪声)
y=Phi*x;                                          %  获得线性测量 
%%  3.  正交匹配追踪法重构信号(本质上是L_1范数最优化问题)
m=2*K;                                            %  算法迭代次数(m>=K)
% Psi=fft(eye(N,N))/sqrt(N);                        %  傅里叶正变换矩阵
for kk=2:N
    for nn=1:N
        dctbasis(kk,nn)=(2/N)^0.5*cos((2*(nn-1)+1)*(kk-1)*pi/2/N);
    end
end
for nn=1:N
    dctbasis(1,nn)=(1/N)^0.5*cos((2*(nn-1)+1)*(1-1)*pi/2/N);
end

Psi=dctbasis;                        %  傅里叶正变换矩阵
T=Phi*Psi';                                       %  恢复矩阵(测量矩阵*正交反变换矩阵)
[r,hat_y,numIts] = newtonromp(K, T,y);
hat_x=real(Psi'*hat_y);                           %  做逆傅里叶变换重构得到时域信号

%%  4.  恢复信号和原始信号对比
error=norm(hat_x-x)^2/norm(x)^2                   %  重构误差
snr=10*log10(1/error)
iter_number=numIts
r=norm(r)

figure(1)
hold on;
plot(hat_x,'k.-')                                 %  重建信号
plot(x,'r')                                       %  原始信号

legend('Recovery','Original')
hold off;