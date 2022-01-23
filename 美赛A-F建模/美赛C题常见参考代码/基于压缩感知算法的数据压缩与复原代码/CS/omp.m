
%Y:orignal signal
%K:the sparse
function [r_n,hat_x,k]=omp(T,N,y,k)
% Phi=randn(M,N);                                 %  测量矩阵(高斯分布白噪声)
% 
%                                            %x是列向量   
% y=Phi*x;                                        %  获得线性测量 
% 
% %    正交匹配追踪法重构信号(本质上是L_1范数最优化问题)
% 
% Psi=fft(eye(N,N))/sqrt(N);                        %  傅里叶正变换矩阵
% T=Phi*Psi';                                       %  恢复矩阵(测量矩阵*正交反变换矩阵)
% 
hat_x=zeros(N,1);                                 %  待重构的谱域(变换域)向量                     
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_n=y;                                            %  残差值

for times=1:k;                                    %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:N;                                  %  恢复矩阵的所有列向量
        product(col)=abs(T(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置
    Aug_t=[Aug_t,T(:,pos)];                       %  矩阵扩充
%     T(:,pos)=zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零）
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;           %  最小二乘,使残差最小
    r_n=y-Aug_t*aug_y;                            %  残差
%     if norm(r_n)<0.001
%         break;
%     end
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
end
hat_x(pos_array)=aug_y;                           %  重构的谱域向量
Iomp=length(pos_array);
% hat_x=real(Psi'*hat_y.');                         %  做逆傅里叶变换重构得到时域信号










