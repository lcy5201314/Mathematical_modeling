%% ����һ����·��·������
% ���룺
% D                 ��������֮��ľ���
% route             һ����·
% �����
% len               ������·����
function len=PathLength(D,route)
[row,col]=size(D);
p=[route route(1)];
i1=p(1:end-1);
i2=p(2:end);
len=sum(D((i1-1)*col+i2));

