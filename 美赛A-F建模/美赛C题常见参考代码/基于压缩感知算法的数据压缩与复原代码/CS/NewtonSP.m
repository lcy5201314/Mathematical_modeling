% %牛顿与sp原子选择结合
% 
% % function x_recovery=Steepest(M,N,x,m)
% clc
% clear
% K=7;      %  稀疏度(做FFT可以看出来)
% N=320;    %  信号长度
% M=100;     %  测量数(M>=K*log(N/K),至少40,但有出错的概率)
% % f1=50;    %  信号频率1
% % f2=100;   %  信号频率2
% % f3=200;   %  信号频率3
% % f4=400;   %  信号频率4
% % fs=800;   %  采样频率
% % ts=1/fs;  %  采样间隔
% % Ts=1:N;   %  采样序列
% % x=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.1*cos(2*pi*f3*Ts*ts)+0.9*cos(2*pi*f4*Ts*ts);  %  完整信号
% % x=x';
% [x,fs,bits]=WAVREAD('M1_1',[8400 8719]);
% A=randn(M,N);
% y=A*x;
% for kk=2:N
%     for nn=1:N
%         dctbasis(kk,nn)=(2/N)^0.5*cos((2*(nn-1)+1)*(kk-1)*pi/2/N);
%     end
% end
% for nn=1:N
%     dctbasis(1,nn)=(1/N)^0.5*cos((2*(nn-1)+1)*(1-1)*pi/2/N);
% end
% 
% B=dctbasis;            
% T=A*B'; %y=A*x=A*B'*s=T*s;    s=B*x;
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
A = randn(M,N);                         %观测矩阵
T=orth(A')';
x=f;
y=T*x;                                          %  获得线性测量 

%初始化
Aug_t=[];
iter_num = 0;
actset_size = k;
active_set = [];
res =y;
m=5;
xsp=zeros(N,1);%恢复出的稀疏信号
% hat_x=zeros(N,1);

% while norm(res)>10^-5
    
    for times=1:20
    [val, idx] = sort(abs(T'*res), 'descend');
    
    candidate_set = union(active_set, idx(1:actset_size));%侯选集合
    a=size(candidate_set,1);
    
    Aug_t=T(:,candidate_set);
    %%%%%%%%%%%%%%%%%%%%%%%%%最速下降
    %     gn=Aug_t'*res;                                   %i*1,梯度
    %
    %     cn=Aug_t*gn;                                     %M*1
    %     d=(res'*cn)/norm(cn).^2;                         %步长
    %     hat_x=hat_x+d*gn  ;
    %     res=res-d*Aug_t*gn;
    
    H=-Aug_t'*Aug_t;                          
    gn=Aug_t'*res;                                
    d=-H^(-1)*gn;                                     
    a=1;
%      xsp( candidate_set)=xsp( candidate_set)+a*d  ;
    
      
   [val idx] = sort(abs(a*d), 'descend');
    new_active_set = candidate_set(idx(1:actset_size));%支撑集合
%     xsp(candidate_set(idx(actset_size:end)))=zeros(actset_size);
    Aug_t=T(:,new_active_set);
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
    active_set= new_active_set;
    H=-Aug_t'*Aug_t;                              %海森矩阵i*i
    gn=Aug_t'*res;                                 %i*1,梯度
    d=-H^(-1)*gn;                                      %步长i*1，d=(rn'*cn)/norm(cn).^2;
    a=1;
        xsp(active_set)=xsp(active_set)+a*d  ;                                              %hat_x(pos_array)=hat_x(pos_array)+d*gn
     res=res-Aug_t* a*d;                                  %a*Aug_t*d
      iter_num = iter_num +1;
end
% xsp(active_set)=xsp(active_set)+a*d  ; 

xrecovery=real(B'*xsp);
iter_num
norm(res)

plot(x,'r');
hold on
plot(xrecovery,'k.-');
hold on;
legend('x','x_recovery');

error=norm(xrecovery-x)^2/norm(x)^2  ;
snr=10*log10(1/error)
hold off


