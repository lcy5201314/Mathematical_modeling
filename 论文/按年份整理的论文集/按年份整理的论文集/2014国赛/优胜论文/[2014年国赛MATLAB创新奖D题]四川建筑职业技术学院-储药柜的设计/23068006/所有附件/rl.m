function [dd,d]=rl(B)%%%%%%B ��ʾrenyiryd�Ľ��,dd �������Ľ��
[m,n]=size(B);

d=tabulate(B(:,2));
d=d(10:end,:);         %%%%%�ۼƹ���������

dd=sortrows(d,3);
