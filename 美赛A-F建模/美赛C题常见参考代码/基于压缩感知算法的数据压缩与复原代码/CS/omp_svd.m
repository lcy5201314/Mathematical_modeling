
%Y:orignal signal
%K:the sparse
function x=omp_svd(Dic,y)
                               
[M,N]=size(Dic);
x=zeros(N,1);                       
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_n=y;                                            %  残差值

for times=1:3                                   %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:N;                                  %  恢复矩阵的所有列向量
        product(col)=abs(Dic(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置
    Aug_t=[Aug_t,Dic(:,pos)];                       %  矩阵扩充
    Dic(:,pos)=zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零）
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;           %  最小二乘,使残差最小
    r_n=y-Aug_t*aug_y;                          %  残差
   
    
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
end
x(pos_array)=aug_y;                           %  重构的稀疏稀疏










