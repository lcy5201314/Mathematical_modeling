
%Y:orignal signal
%K:the sparse
function x=omp_svd(Dic,y)
                               
[M,N]=size(Dic);
x=zeros(N,1);                       
Aug_t=[];                                         %  ��������(��ʼֵΪ�վ���)
r_n=y;                                            %  �в�ֵ

for times=1:3                                   %  ��������(�������������,�õ�������ΪK)
    for col=1:N;                                  %  �ָ����������������
        product(col)=abs(Dic(:,col)'*r_n);          %  �ָ�������������Ͳв��ͶӰϵ��(�ڻ�ֵ) 
    end
    [val,pos]=max(product);                       %  ���ͶӰϵ����Ӧ��λ��
    Aug_t=[Aug_t,Dic(:,pos)];                       %  ��������
    Dic(:,pos)=zeros(M,1);                          %  ѡ�е������㣨ʵ����Ӧ��ȥ����Ϊ�˼��Ұ������㣩
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;           %  ��С����,ʹ�в���С
    r_n=y-Aug_t*aug_y;                          %  �в�
   
    
    pos_array(times)=pos;                         %  ��¼���ͶӰϵ����λ��
end
x(pos_array)=aug_y;                           %  �ع���ϡ��ϡ��










