
%Y:orignal signal
%K:the sparse
function [r_n,hat_x,k]=omp(T,N,y,k)
% Phi=randn(M,N);                                 %  ��������(��˹�ֲ�������)
% 
%                                            %x��������   
% y=Phi*x;                                        %  ������Բ��� 
% 
% %    ����ƥ��׷�ٷ��ع��ź�(��������L_1�������Ż�����)
% 
% Psi=fft(eye(N,N))/sqrt(N);                        %  ����Ҷ���任����
% T=Phi*Psi';                                       %  �ָ�����(��������*�������任����)
% 
hat_x=zeros(N,1);                                 %  ���ع�������(�任��)����                     
Aug_t=[];                                         %  ��������(��ʼֵΪ�վ���)
r_n=y;                                            %  �в�ֵ

for times=1:k;                                    %  ��������(�������������,�õ�������ΪK)
    for col=1:N;                                  %  �ָ����������������
        product(col)=abs(T(:,col)'*r_n);          %  �ָ�������������Ͳв��ͶӰϵ��(�ڻ�ֵ) 
    end
    [val,pos]=max(product);                       %  ���ͶӰϵ����Ӧ��λ��
    Aug_t=[Aug_t,T(:,pos)];                       %  ��������
%     T(:,pos)=zeros(M,1);                          %  ѡ�е������㣨ʵ����Ӧ��ȥ����Ϊ�˼��Ұ������㣩
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;           %  ��С����,ʹ�в���С
    r_n=y-Aug_t*aug_y;                            %  �в�
%     if norm(r_n)<0.001
%         break;
%     end
    pos_array(times)=pos;                         %  ��¼���ͶӰϵ����λ��
end
hat_x(pos_array)=aug_y;                           %  �ع�����������
Iomp=length(pos_array);
% hat_x=real(Psi'*hat_y.');                         %  ���渵��Ҷ�任�ع��õ�ʱ���ź�










