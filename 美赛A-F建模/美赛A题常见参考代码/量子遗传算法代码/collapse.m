function binary=collapse(chrom)
%% ����Ⱥʵʩһ�β��� �õ������Ʊ���
% ����chrom ��Ϊ���ӱ��ر���
% ���binary�������Ʊ���
[M,N]=size(chrom);  %�õ���Ⱥ��С �ͱ��볤��
M=M/2;  % ��Ⱥ��С
binary=zeros(M,N);  %�����Ʊ����С��ʼ��
for i=1:M
    for j=1:N
        pick=rand;  %������0,1�������
        if pick>(chrom(2.*i-1,j)^2)    % ��������ڦ���ƽ��
            binary(i,j)=1;
        else
            binary(i,j)=0;
        end
    end
end