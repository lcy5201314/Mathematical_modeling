% %牛顿与sp原子选择结合
% 
% % function x_recovery=Steepest(M,N,x,m)
% clc
% clear
% k=16;      %  稀疏度(做FFT可以看出来)
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

%初始化
tic
clc
clear
% N is the original signal length
N = 2^10;
% k is number of observations to make
M = 450;

% number of spikes to put down
% n_spikes = floor(.01*n);
k= 100;


% random +/- 1 signal
f = zeros(N,1);
q = randperm(N);
f(q(1:k)) = randn(k,1);
disp('Creating measurement matrix...');
A = randn(M,N);                         %观测矩阵
T=orth(A')';
x=f;

y=T*x;
numIts = 0;

iter_num = 0;
actset_size =k;
active_set = [];
New_active_set=[];
res =y;

xr=zeros(N,1);
% hat_x=zeros(N,1);


    
   
%Run ROMP
I=[];
Jnew=[];
    
while length(I)-1 <2*k && norm(res) > 10^(-6)
      
        %Find J, the biggest n coordinates of u
        u = T' * res;
        absu = abs(u);
        [b, ix] = sort(absu, 'descend');
        J = ix(1:k);
        Jvals = b(1:k);
        %Find J0, the set of comparable coordinates with maximal energy
        w=1;
        best = -1;
        J0 = zeros(1);
        while w <= k
            first = Jvals(w);
            firstw = w;
            energy = 0;
            while ( w <= k ) && ( Jvals(w) >= 1/2 * first )
                energy = energy + Jvals(w)^2;
                w = w+1;
            end
            if energy > best
                best = energy;
                J0 = J(firstw:w-1);
            end
        end
        %Add J0 to the index set I
        candidate_set=ix(1:actset_size);
         Aug_t=T(:,J0);
          hat_x=pinv(Aug_t)*y;
        [val idx] = sort(abs(hat_x), 'descend');
         Jnew = candidate_set(idx(1:length(J0)));
         
         
        
        I(length(I)+1: length(I)+length(Jnew)) = Jnew;
        I=unique(I);
       
        
        
       
       
        Aug_t=T(:,I);
        res= y - Aug_t*pinv(Aug_t) * y;
       
        iter_num=iter_num+1;
end
 
%  for times=1:5
%      if mod(times,2)==1
%     candidate_set = union(active_set, I);
%      else 
%          u = T' * res;
%         absu = abs(u);
%         [b, ix] = sort(absu, 'descend');
%         candidate_set=union(active_set,ix(1:actset_size));
%      end
%        
%         
%         Aug_t=T(:,candidate_set);
%         
%         hat_x=pinv(Aug_t)*y;
%         
%         
%         [val idx] = sort(abs(hat_x), 'descend');
%         new_active_set = candidate_set(idx(1:actset_size));%回溯（从侯选集合选支撑集合）
%         active_set= new_active_set;
%         Aug_t=T(:,new_active_set);
%         new_res = y-Aug_t*pinv(Aug_t)*y;
%          res = new_res;

vSmall = Aug_t\ y;
vOut = zeros(N, 1);
vOut(I)=vSmall;
I_length=length(I)
% xr=real(B'*xr);

iter_num
rn=norm(res)

plot(x,'r');
hold on
% plot(xr,'k.-');
plot(vOut,'k.-')
hold on;
legend('x','xr');

error=norm(vOut-x)^2/norm(x)^2  ;
snr=10*log10(1/error)
hold off
toc

