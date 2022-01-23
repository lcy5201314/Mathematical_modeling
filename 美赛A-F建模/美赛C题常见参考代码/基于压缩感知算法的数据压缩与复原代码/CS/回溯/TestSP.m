clc;
clear;

%  1. ʱ������ź�����
K=30;       %  ϡ���(��FFT���Կ�����)
N=320;      %  �źų���
M=100;      %  ������(M>=K*log(N/K),����40,���г���ĸ���)
m=18;
% f1=50;    %  �ź�Ƶ��1
% f2=100;   %  �ź�Ƶ��2
% f3=200;   %  �ź�Ƶ��3
% f4=400;   %  �ź�Ƶ��4
% fs=800;   %  ����Ƶ��
% ts=1/fs;  %  �������
% Ts=1:N;   %  ��������
% x=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.1*cos(2*pi*f3*Ts*ts)+0.9*cos(2*pi*f4*Ts*ts);  %  �����ź�
% x=x';
[x,fs,bits]=wavread('M1_1',[8400 8719]);
%  2.  ʱ���ź�ѹ������
Psi=fft(eye(N,N))/sqrt(N);                        %  ����Ҷ���任����
   
Phi=randn(M,N);  
T=Phi*Psi';                                       %  ��������(��˹�ֲ�������)
y=Phi*x;                                          %  ������Բ��� 
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


