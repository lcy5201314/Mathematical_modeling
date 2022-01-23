%  1-D信号压缩传感的实现(正交匹配追踪法Orthogonal Matching Pursuit)
%  测量数M>=K*log(N/K),K是稀疏度,N信号长度,可以近乎完全重构
%  编程人--香港大学电子工程系 沙威  Email: wsha@eee.hku.hk
%  编程时间：2008年11月18日
%  文档下载: http://www.eee.hku.hk/~wsha/Freecode/freecode.htm 
%  参考文献：Joel A. Tropp and Anna C. Gilbert 
%  Signal Recovery From Random Measurements Via Orthogonal Matching
%  Pursuit，IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 53, NO. 12,
%  DECEMBER 2007.
tic
clc;
clear;

%%
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
%%

Phi=randn(M,N);                                   %  测量矩阵(高斯分布白噪声)
y=Phi*x;                                          %  获得线性测量 y=phi*x,s=psi*x =>y=phi*psi'*s=T*s
s=zeros(N,1);                                     %  待重构的谱域(变换域)向量                     
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_n=y;  

%%
m=2*K;                                            %  算法迭代次数(m>=K)
Psi=fft(eye(N,N))/sqrt(N);                        %  傅里叶正变换矩阵
T=Phi*Psi';                                       %  恢复矩阵(测量矩阵*正交反变换矩阵)
for times=1:12                                   %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:N                                  %  恢复矩阵的所有列向量
        product(col)=abs(T(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);    
    pos_array(times)=pos;                       %  最大投影系数对应的位置
    valnew=val;
    Aug_t=[Aug_t,T(:,pos)];                       %  矩阵扩充
                                                 
    s(pos)=s(pos)+valnew;
    r_n=r_n-valnew*T(:,pos)./norm(T(:,pos));                      %  残差
                             %  选中的列置零（实质上应该去掉，为了简单我把它置零）
                           %  纪录最大投影系数的位置
end


          
hat_x=real(Psi'*s);   %  做逆傅里叶变换重构得到时域信号
plot(hat_x,'k.-');
hold on;
plot(x,'r')
legend('Recovery','Original')
toc