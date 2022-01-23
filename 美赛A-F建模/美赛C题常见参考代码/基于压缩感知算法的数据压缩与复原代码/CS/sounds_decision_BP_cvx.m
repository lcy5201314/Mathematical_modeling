% 用convex工具箱来作BP算法
clear all;
close all;
clc;
cvx_quiet(true);%可以先安装cvx_setup

%E1=[];
E2=[];
%E3=[];
%for snr=5:5:30              %不同信噪比下的信号恢复
N=320;
% X=wavread('F:/academic/code/sounds/SYL1.wav');
X=wavread('F:/academic/code/sounds/SYL1.wav',[8000,9000]);
x=enframe(X,N);
x=x';
    [N,L]=size(x);
    M=150;
    Psi=fft(eye(N,N))/sqrt(N);                        %  傅里叶正变换矩阵
    hat_x=zeros(N,L);
for k=1:L
    y=Psi*x(:,k);                                          %  K稀疏信号
    y_real=real(y);
    y_imag=imag(y);
    
    Phi = randn(M,N);% Measurement matrix
    T=Phi*Psi';
    u_real = T*y_real;      %Generate measurements
    %snr=20;         %dB
    %delta=10^(-snr/20);
    %u=u+delta*randn(size(u));       %加噪声
    %-------BP----
    cvx_begin
    variable hat_y_real(N,1);
    minimize (norm(hat_y_real,1));
    subject to
    u_real==T*hat_y_real;
    cvx_end
    
    u_imag = T*y_imag;      %Generate measurements
    %snr=20;         %dB
    %delta=10^(-snr/20);
    %u=u+delta*randn(size(u));       %加噪声
    %-------BP----
    cvx_begin
    variable hat_y_imag(N,1);
    minimize (norm(hat_y_imag,1));
    subject to
    u_imag==T*hat_y_imag;
    cvx_end

    hat_y=hat_y_real+i*hat_y_imag;
    hat_x(:,k)=real(Psi'*hat_y);
%     e2=norm(hat_y-y,2);
%     E2=[E2,e2];
end
hat_X=overlapadd(hat_x');
% snr=5:5:30;
% figure(1)
% plot(snr,E2,'r-*');hold on;
% grid on;
% xlabel('SNR(dB)');
% ylabel('MSE');

% figure(2)
% plot(abs(hat_y),'r');
% hold on
% plot(abs(y),'b');
% hold off


%hat_x=fliplr(hat_x);                           %信号翻转
figure(3)
hold on;
plot(hat_X,'k.-')                                 %  重建信号
plot(X,'r')                                       %  原始信号
legend('Recovery','Original')
hold off
% norm(hat_x-x)/norm(x)                           %  重构误差