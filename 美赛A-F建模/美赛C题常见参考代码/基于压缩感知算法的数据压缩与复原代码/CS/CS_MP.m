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
clc;
clear;

%%
K=7;      %  ϡ���(��FFT���Կ�����)
N=256;    %  �źų���
M=64;     %  ������(M>=K*log(N/K),����40,���г���ĸ���)
f1=50;    %  �ź�Ƶ��1
f2=100;   %  �ź�Ƶ��2
f3=200;   %  �ź�Ƶ��3
f4=400;   %  �ź�Ƶ��4
fs=800;   %  ����Ƶ��
ts=1/fs;  %  �������
Ts=1:N;   %  ��������
x=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.1*cos(2*pi*f3*Ts*ts)+0.9*cos(2*pi*f4*Ts*ts);  %  �����ź�
x=x';
%%

Phi=randn(M,N);                                   %  ��������(��˹�ֲ�������)
y=Phi*x;                                          %  ������Բ��� y=phi*x,s=psi*x =>y=phi*psi'*s=T*s
s=zeros(N,1);                                     %  ���ع�������(�任��)����                     
Aug_t=[];                                         %  ��������(��ʼֵΪ�վ���)
r_n=y;  

%%
m=2*K;                                            %  �㷨��������(m>=K)
Psi=fft(eye(N,N))/sqrt(N);                        %  ����Ҷ���任����
T=Phi*Psi';                                       %  �ָ�����(��������*�������任����)
for times=1:12                                   %  ��������(�������������,�õ�������ΪK)
    for col=1:N                                  %  �ָ����������������
        product(col)=abs(T(:,col)'*r_n);          %  �ָ�������������Ͳв��ͶӰϵ��(�ڻ�ֵ) 
    end
    [val,pos]=max(product);    
    pos_array(times)=pos;                       %  ���ͶӰϵ����Ӧ��λ��
    valnew=val;
    Aug_t=[Aug_t,T(:,pos)];                       %  ��������
                                                 
    s(pos)=s(pos)+valnew;
    r_n=r_n-valnew*T(:,pos)./norm(T(:,pos));                      %  �в�
                             %  ѡ�е������㣨ʵ����Ӧ��ȥ����Ϊ�˼��Ұ������㣩
                           %  ��¼���ͶӰϵ����λ��
end


          
hat_x=real(Psi'*s);   %  ���渵��Ҷ�任�ع��õ�ʱ���ź�
plot(hat_x,'k.-');
hold on;
plot(x,'r')
legend('Recovery','Original')
toc