tic
clc
clear
% N is the original signal length
N = 2^10;
% k is number of observations to make
M = 450;
% number of spikes to put down
% n_spikes = floor(.01*n);
k= 120;
% random +/- 1 signal
f = zeros(N,1);
q = randperm(N);
f(q(1:k)) = randn(k,1);
disp('Creating measurement matrix...');
A = randn(M,N);                         %观测矩阵
T=orth(A')';
x=f;
y=T*x;
figure(1)
plot(x,'r');
hold on;
Aug_t=[];
rn=y;
hat_x=zeros(N,1);
%%%%%%%%%%%%%%%runnewton
 for times=1:120
    for col=1:N
        inner(col)=abs(T(:,col)'*rn);
    end
    
    [val,pos]=max(inner);
    Aug_t=[Aug_t,T(:,pos)];  %M*i
   l=size(Aug_t,2);
    pos_array(times)=pos; %i*1
    

     H=-Aug_t'*Aug_t;                              %海森矩阵i*i
    gn=Aug_t'*rn;                                 %i*1,梯度
     d=-H^(-1)*gn;                                      %步长i*1，d=(rn'*cn)/norm(cn).^2;
%     a=-(gn'*d)/(d'*H*d);
     hat_x(pos_array)=hat_x(pos_array)+d  ;         %hat_x(pos_array)=hat_x(pos_array)+d*gn  
   
  
   rn=rn-Aug_t*d;                                  %a*Aug_t*d
     T(:,pos)=zeros(M,1); 
% if norm(rn)<10^-3
%      break;
%  end
% times=times+1;
   
end
norm(rn)
iter_num=times
error=norm(hat_x-x)^2/norm(x)^2  ;
snr=10*log10(1/error)

disp('draw the recovery signal ...')
plot(hat_x,'k.-');
hold off
toc