% ����ͼ�Ĺ���������ڽӾ���ת��
function W=mattransf(G,f)
% f=0, �ڽӾ���ת��Ϊ��������
% f=1����������ת��Ϊ�ڽӾ���

if f==0 % �ڽӾ���ת��Ϊ��������
    m=sum(sum(G));
    n=size(G,1);
    W=zeros(n,m);
    k=1;
    for i=1:n
        for j=1:n
            if G(i,j)~=0 % ��i�����ıߣ�����ߵ�ʼ��
                W(i,k)=1; % ��������ʼ��ֵΪ1
                W(j,k)=-1; % ���������յ�ֵΪ-1
                k=k+1;
            end
        end
    end
elseif f==1 % ��������ת��Ϊ�ڽӾ���
    m=size(G,2);
    n=size(G,1);
    W=zeros(n,n);
    for i=1:m
        a=find(G(:,i)~=0); % ����ߵ���������
        if G(a(1),i)==1
            W(a(1),a(2))=1; % �������a(1)ָ��a(2)
        else
            W(a(2),a(1))=1; % �������a(2)ָ��a(1)
        end
    end
else
    fprint('please input the right value of f');
end
W;
