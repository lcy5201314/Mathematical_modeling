clc;clear

%%  1. 时域测试信号生成
K=7;      %  稀疏度(做FFT可以看出来)
N=256;    %  信号长度
M=64;     %  测量数(M>=K*log(N/K),至少40,但有出错的概率)
f1=50;    %  信号频率1
f2=100;   %  信号频率2
f3=200;   %  信号频率3
f4=400;   %  信号频率4
fs=800;   %  采样频率
ts=1/fs;  %  采样间隔
Ts=1:N;   %  采样序列
x1=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.1*cos(2*pi*f3*Ts*ts)+0.9*cos(2*pi*f4*Ts*ts);  %  完整信号
tau=ones(64,1);

%%  2.  时域信号压缩传感
Phi=randn(M,N);                                   %  测量矩阵(高斯分布白噪声)
y=Phi*x1.';%  获得线性测量
A=@(x1) Phi*x1;
At=@(x1) Phi'*x1;
x=;
figure(1);
plot(x1,'k.')
hold on;
plot(y,'r.')
%[x,x_debias,objective,times,debias_start,mses]=GPSR_Basic(y,A,tau,'Debias',1,...
 %        'AT',At,... 
  %       'Initialization',0,...
   % 	 'StopCriterion',3,...     'ToleranceA',0.01,...
    %     'ToleranceD',0.0001);
figure(2);
plot(x1,'k-')
hold on;
plot(x,'r.-')
