%%
%���ٶ��½�����spԭ��ѡ����

% function x_recovery=Steepest(M,N,x,m)
clc
clear
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
A=randn(M,N);
y=A*x;
B=fft(eye(N,N))/sqrt(N); %y=A*x=A*B'*s=T*s;    s=B*x;
T=A*B';
%%

%��ʼ��
Aug_t=[];
iter_num = 0;
actset_size = K;
active_set = [];
res =y;
m=5;
xsp=zeros(N,1);%�ָ�����ϡ���ź�
% hat_x=zeros(N,1);

% while norm(res)>10^-5
    
    for times=1:40
    [val, idx] = sort(abs(T'*res), 'descend');
    
    candidate_set = union(active_set, idx(1:actset_size));%��ѡ����
    a=size(candidate_set,1);
    
    Aug_t=T(:,candidate_set);
    %%%%%%%%%%%%%%%%%%%%%%%%%�����½�
        gn=Aug_t'*res;                                   %i*1,�ݶ�
    
        cn=Aug_t*gn;                                     %M*1
        d=(res'*cn)/norm(cn).^2;                         %����
      
   
    

%      xsp( candidate_set)=xsp( candidate_set)+a*d  ;
    
      
   [val idx] = sort(abs(d*gn), 'descend');
    new_active_set = candidate_set(idx(1:actset_size));%֧�ż���
%     xsp(candidate_set(idx(actset_size:end)))=zeros(actset_size);
    Aug_t=T(:,new_active_set);
    active_set= new_active_set;
%     %      new_res=y-Aug_t*pinv(Aug_t)*y;%��С���˻ָ�
%     %
%     %       res = new_res;
%     %      active_set= new_active_set;
%     %
%     %     iter_num = iter_num +1;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ʹ�������½���
%     %
%     %     gn=Aug_t'*res;                                   %i*1,�ݶ�
%     %     cn=Aug_t*gn;                                     %M*1
%     %     d=(res'*cn)/norm(cn).^2;                         %����
%     %
%     %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %     res=res-d*Aug_t*gn;
%     %
%     %     active_set= new_active_set;
%     %      xsp(active_set) = xsp(active_set)+d*gn;
%     %    
%     %     T(:,active_set)=zeros(M,size(active_set,1));
          gn=Aug_t'*res;                                   %i*1,�ݶ�Ҳ����Ѱ�ŷ���
    
        cn=Aug_t*gn;                                     %M*1
        d=(res'*cn)/norm(cn).^2;                         %����
      
   
        xsp(active_set)=xsp(active_set)+d*gn  ;                                              %hat_x(pos_array)=hat_x(pos_array)+d*gn
     res=res-d*Aug_t* gn;                                  %a*Aug_t*d
      iter_num = iter_num +1;
end
% xsp(active_set)=xsp(active_set)+a*d  ; 

xrecovery=real(B'*xsp);
iter_num
norm(res)
plot(xrecovery,'k.-');
hold on;
plot(x,'r');
legend('x','x_recovery')
error=norm(xrecovery-x)^2/norm(x)^2  ;
snr=10*log10(1/error)
hold off


