%% �����˹���Ⱥ��ǰλ�õ�ʳ��Ũ��
%����neighbork                ��Ⱥ����
%����D��                      �������
%���Y��                      �˹���Ⱥ����·�����ܾ���
function Y=AF_foodconsistence(neighbork,D)

neigh_Num=size(neighbork,1);                %��Ⱥ��Ŀ
Y=zeros(neigh_Num,1);
for i=1:neigh_Num
    Y(i)=PathLength(D,neighbork(i,:));
end


end

