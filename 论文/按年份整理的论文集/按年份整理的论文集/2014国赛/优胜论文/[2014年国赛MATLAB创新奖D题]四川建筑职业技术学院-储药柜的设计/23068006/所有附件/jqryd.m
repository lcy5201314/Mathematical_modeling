function [Q,P]=jqryd(A,B)%%%%��Ȩ����ȼ��㣬A��ʾrenyiryd��������Ľ����B��ʾrl��������������dd,P��ʾ��Ȩ��������,Q �ڶ����Ǽ�Ȩ�����
[m,n]=size(B);
for i=1:m
    b=B(i,1);
    r=min(A(A(:,2)==b,4));
    Q(i,2)=B(i,3)*r;
    Q(i,1)=b;
    Q(i,3)=r;
    Q(i,4)=B(i,3);
end
   Q=sortrows(Q,2);
    P=sum(Q(:,2));