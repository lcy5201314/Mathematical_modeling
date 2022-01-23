clc;
clear;

%%  1. ʱ������ź�����
K=10;      %  ϡ���(��FFT���Կ�����)
N=320;    %  �źų���
M=100;     %  ������(M>=K*log(N/K),����40,���г���ĸ���)
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
%%  2.  ʱ���ź�ѹ������
Phi=randn(M,N);                                   %  ��������(��˹�ֲ�������)
y=Phi*x;                                          %  ������Բ��� 
%%  3.  ����ƥ��׷�ٷ��ع��ź�(��������L_1�������Ż�����)
m=2*K;                                            %  �㷨��������(m>=K)
% Psi=fft(eye(N,N))/sqrt(N);                        %  ����Ҷ���任����
for kk=2:N
    for nn=1:N
        dctbasis(kk,nn)=(2/N)^0.5*cos((2*(nn-1)+1)*(kk-1)*pi/2/N);
    end
end
for nn=1:N
    dctbasis(1,nn)=(1/N)^0.5*cos((2*(nn-1)+1)*(1-1)*pi/2/N);
end

Psi=dctbasis;                        %  ����Ҷ���任����
T=Phi*Psi';                                       %  �ָ�����(��������*�������任����)
[r,hat_y,numIts] = newtonromp(K, T,y);
hat_x=real(Psi'*hat_y);                           %  ���渵��Ҷ�任�ع��õ�ʱ���ź�

%%  4.  �ָ��źź�ԭʼ�źŶԱ�
error=norm(hat_x-x)^2/norm(x)^2                   %  �ع����
snr=10*log10(1/error)
iter_number=numIts
r=norm(r)

figure(1)
hold on;
plot(hat_x,'k.-')                                 %  �ؽ��ź�
plot(x,'r')                                       %  ԭʼ�ź�

legend('Recovery','Original')
hold off;