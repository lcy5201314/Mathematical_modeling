function C=chanel_cap(P,e)
[r,s]=size(P);
n=0;
C=0;
C_0=0;
C_1=0;
X=ones(1,r)/r;
A=zeros(1,r);
B=zeros(r,s);
while(1)
    n=n+1;
    for i=1:r
        for j=1:s

          
         if P(i,j)==0 
                B(i,j)=0;
         else
                  B(i,j)=log(P(i,j)/(X*P(:,j)));
            end


        end  
        A(1,i)=exp(P(i,:)*B(i,:)');       
    end    
    C_0=log2(X*A');
    C_1=log2(max(A));
    if (abs(C_0-C_1)<e)  %满足迭代终止条件停止迭代 
        C=C_0;        
        fprintf('迭代次数: n=%d\n',n)
        fprintf('信道容量: C=%f比特/符号',C)
        break;  %满足后输出结果并退出
     else
         X=(X.*A)/(X*A');
         continue;        
     end
end
