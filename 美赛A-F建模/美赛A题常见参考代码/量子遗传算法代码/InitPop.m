function chrom=InitPop(M,N)
%% ��ʼ����Ⱥ-���ӱ��ر���
% M:Ϊ��Ⱥ��С��2��(���ͦ�)
% N:Ϊ���ӱ��ر��볤��
for i=1:M
    for j=1:N
        chrom(i,j)=1/sqrt(2);
    end
end