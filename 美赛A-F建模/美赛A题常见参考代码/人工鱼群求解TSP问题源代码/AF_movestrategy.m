%% �ƶ�����
%�Ƚ���׷β��Ϊ,���û�н����ٽ�����ʳ��Ϊ,
%�����û�н�������о�Ⱥ��Ϊ,�����Ȼû�н�����
%��������ƶ���Ϊ.
%����X��               ��Ⱥ����
%����i��               ��i���˹���
%����D��               �������
%����Visual��          ��֪����
%����deta��            ӵ��������
%����trynumber��       �����̽����
%���Xinext��          ���ҵ���·��
%���flag��            ����Ƿ��ҵ����õ�·����flag=0��ʾ�½�û�иĽ���flag=1��ʾ�½��иĽ�
function [Xinext,flag]= AF_movestrategy(X,i,D,Visual,deta,trynumber)
Xi=X(i,:);                                                              %��i���˹���
Yi=PathLength(D,Xi);                                                    %·��Xi���ܾ���
flag=0;
flag_prey=1;
flag_swarm=1;

%% ׷β��Ϊ
[Xinext,flag_follow]=AF_follow(X,i,D,Visual,deta,trynumber);

%% ���û�н����ٽ�����ʳ��Ϊ
if flag_follow==0
    [Xinext,flag_prey]=AF_prey(X,i,D,trynumber,Visual);
end

%% �����û�н�������о�Ⱥ��Ϊ
if flag_prey==0
    [Xinext,flag_swarm]=AF_swarm(X,i,D,Visual,deta,trynumber);
end

%% �����Ȼû�н����ͽ�������ƶ���Ϊ
if flag_swarm==0
    Xinext = AF_randmove(Xi);
end

%% ����ж��½�Xinext�Ƿ��иĽ�
Yinext=PathLength(D,Xinext);                                               %·��Xinext���ܾ���
if Yinext<Yi
    flag=1;
end

end

