% clc;
% clear;
% close
% %  1. ʱ������ź�����
% K=16;       %  ϡ���(��FFT���Կ�����)
% N=320;      %  �źų���
% M=100;      %  ������(M>=K*log(N/K),����40,���г���ĸ���)
% m=2;%��������
% % f1=50;    %  �ź�Ƶ��1
% % f2=100;   %  �ź�Ƶ��2
% % f3=200;   %  �ź�Ƶ��3
% % f4=400;   %  �ź�Ƶ��4
% % fs=800;   %  ����Ƶ��
% % ts=1/fs;  %  �������
% % Ts=1:N;   %  ��������
% % x=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.1*cos(2*pi*f3*Ts*ts)+0.9*cos(2*pi*f4*Ts*ts);  %  �����ź�
% % x=x';
% [x,fs,bits]=wavread('M1_1',[8400 8719]);
% %  2.  ʱ���ź�ѹ������
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
tic%  ��������(��˹�ֲ�������)
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
A = randn(M,N);                         %�۲����
T=orth(A')';
x=f;
y=T*x;                                          %  ������Բ��� 
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

