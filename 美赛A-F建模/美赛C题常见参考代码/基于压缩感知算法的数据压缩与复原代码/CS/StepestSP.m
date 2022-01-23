%%
%最速度下降法与sp原子选择结合

% function x_recovery=Steepest(M,N,x,m)
clc
clear
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
x=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.1*cos(2*pi*f3*Ts*ts)+0.9*cos(2*pi*f4*Ts*ts);  %  完整信号
x=x';
A=randn(M,N);
y=A*x;
B=fft(eye(N,N))/sqrt(N); %y=A*x=A*B'*s=T*s;    s=B*x;
T=A*B';
%%

%初始化
Aug_t=[];
iter_num = 0;
actset_size = K;
active_set = [];
res =y;
m=5;
xsp=zeros(N,1);%恢复出的稀疏信号
% hat_x=zeros(N,1);

% while norm(res)>10^-5
    
    for times=1:40
    [val, idx] = sort(abs(T'*res), 'descend');
    
    candidate_set = union(active_set, idx(1:actset_size));%侯选集合
    a=size(candidate_set,1);
    
    Aug_t=T(:,candidate_set);
    %%%%%%%%%%%%%%%%%%%%%%%%%最速下降
        gn=Aug_t'*res;                                   %i*1,梯度
    
        cn=Aug_t*gn;                                     %M*1
        d=(res'*cn)/norm(cn).^2;                         %步长
      
   
    

%      xsp( candidate_set)=xsp( candidate_set)+a*d  ;
    
      
   [val idx] = sort(abs(d*gn), 'descend');
    new_active_set = candidate_set(idx(1:actset_size));%支撑集合
%     xsp(candidate_set(idx(actset_size:end)))=zeros(actset_size);
    Aug_t=T(:,new_active_set);
    active_set= new_active_set;
%     %      new_res=y-Aug_t*pinv(Aug_t)*y;%最小二乘恢复
%     %
%     %       res = new_res;
%     %      active_set= new_active_set;
%     %
%     %     iter_num = iter_num +1;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%使用最速下降法
%     %
%     %     gn=Aug_t'*res;                                   %i*1,梯度
%     %     cn=Aug_t*gn;                                     %M*1
%     %     d=(res'*cn)/norm(cn).^2;                         %步长
%     %
%     %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %     res=res-d*Aug_t*gn;
%     %
%     %     active_set= new_active_set;
%     %      xsp(active_set) = xsp(active_set)+d*gn;
%     %    
%     %     T(:,active_set)=zeros(M,size(active_set,1));
          gn=Aug_t'*res;                                   %i*1,梯度也就是寻优方向
    
        cn=Aug_t*gn;                                     %M*1
        d=(res'*cn)/norm(cn).^2;                         %步长
      
   
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


