function W=incandadf(G,f)
% ����������ڽӾ����ת��
% G ͼ����Ӧ����
% f=0, �ڽӾ���ת��Ϊ��������
% f=1����������ת��Ϊ�ڽӾ���
% W ת�����

if f==0 % �ڽӾ���ת��Ϊ��������
    m=sum(sum(G))/2; % ����ͼ�ı���
    n=size(G,1);
    W=zeros(n,m);
    k=1;
    for i=1:n
        for j=i:n
            if G(i,j)~=0
                W(i,k)=1; % ���ߵ�ʼ�㸳ֵΪ1
                W(j,k)=1; % ���ߵ��յ㸳ֵΪ1
                k=k+1;
            end
        end
    end
elseif f==1 % ��������ת��Ϊ�ڽӾ���
    m=size(G,2);
    n=size(G,1);
    W=zeros(n,n);
    for i=1:m
        a=find(G(:,i)~=0);
        W(a(1),a(2))=1;  %  ���ڱ����ڽӾ���Ķ�ӦֵΪ1
        W(a(2),a(1))=1;
    end
else
    fprint('please input the right value of f');
end
W;
% G=[0 1 1 1;1 0 1 1;1 1 0 1;1 1 1 0];
% G=[1 1 1 0 0 0;1 0 0 1 1 0;0 1 0 1 0 1;0 0 1 0 1 1];