function X=bin2decFun(x,lenchrom,bound)
%% ������ת����ʮ����
% ����      x�������Ʊ���
%    lenchrom���������Ķ�����λ��
%       bound���������ķ�Χ
% ���      X��ʮ������
M=length(lenchrom);
n=1;
X=zeros(1,M);
for i=1:M
    for j=lenchrom(i)-1:-1:0
        X(i)=X(i)+x(n).*2.^j;
        n=n+1;
    end
end
X=bound(:,1)'+X./(2.^lenchrom-1).*(bound(:,2)-bound(:,1))'; 