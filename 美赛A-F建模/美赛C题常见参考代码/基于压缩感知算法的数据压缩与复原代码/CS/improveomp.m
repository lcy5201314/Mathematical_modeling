%  1-D�ź�ѹ�����е�ʵ��(����ƥ��׷�ٷ�Orthogonal Matching Pursuit)
%  ������M>=K*log(N/K),K��ϡ���,N�źų���,���Խ�����ȫ�ع�
%  �����--��۴�ѧ���ӹ���ϵ ɳ��  Email: wsha@eee.hku.hk
%  ���ʱ�䣺2008��11��18��
%  �ĵ�����: http://www.eee.hku.hk/~wsha/Freecode/freecode.htm 
%  �ο����ף�Joel A. Tropp and Anna C. Gilbert 
%  Signal Recovery From Random Measurements Via Orthogonal Matching
%  Pursuit��IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 53, NO. 12,
%  DECEMBER 2007.
tic
clear
clc
close
[x,fs,bits]=WAVREAD('M1_1',[8400 8719]);
%%  1. ʱ������ź�����
K=7;      %  ϡ���(��FFT���Կ�����)
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

%%  2.  ʱ���ź�ѹ������
Phi=randn(M,N);                                   %  ��������(��˹�ֲ�������)
s=Phi*x;                                        %  ������Բ��� 

%%  3.  ����ƥ��׷�ٷ��ع��ź�(��������L_1�������Ż�����)
m=2*K;                                            %  �㷨��������(m>=K)
Psi=fft(eye(N,N))/sqrt(N);                        %  ����Ҷ���任����
T=Phi*Psi';                                       %  �ָ�����(��������*�������任����)

hat_y=zeros(1,N);                                 %  ���ع�������(�任��)����                     
Aug_t=[];                                         %  ��������(��ʼֵΪ�վ���)
r_n=s;                                            %  �в�ֵ
lamda=1:0.1:20;
    lamda1=lamda(:,1:M);
% lamda1=randn(1,M);
for times=1:m;                                    %  ��������(�������������,�õ�������ΪK)
   
    [rp,index]=sort(r_n,'descend');
    rp1=rp.*lamda1';
    for col=1:M
        r_n(index(col))=rp1(col);
    end
    for col=1:N;                                  %  �ָ����������������
        product(col)=abs(T(:,col)'*r_n);          %  �ָ�������������Ͳв��ͶӰϵ��(�ڻ�ֵ) 
    end
    
    [val,pos]=max(product);                       %  ���ͶӰϵ����Ӧ��λ��
    Aug_t=[Aug_t,T(:,pos)];                       %  ��������
    T(:,pos)=zeros(M,1);                          %  ѡ�е������㣨ʵ����Ӧ��ȥ����Ϊ�˼��Ұ������㣩
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s; %  ��С����,ʹ�в���С
    r_n=s-Aug_t*aug_y;                            %  �в�
    pos_array(times)=pos;                         %  ��¼���ͶӰϵ����λ��
end
hat_y(pos_array)=aug_y;                           %  �ع�����������
hat_x=real(Psi'*hat_y.');                         %  ���渵��Ҷ�任�ع��õ�ʱ���ź�

%%  4.  �ָ��źź�ԭʼ�źŶԱ�
figure(2);
hold on;
plot(hat_x,'k.-')                                 %  �ؽ��ź�
plot(x,'r')                                       %  ԭʼ�ź�
legend('Recovery','Original')
norm(hat_x-x)^2/norm(x)^2                         %  �ع����
toc