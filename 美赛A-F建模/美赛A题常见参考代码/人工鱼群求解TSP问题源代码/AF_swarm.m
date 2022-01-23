%% ��Ⱥ��Ϊ
%����X��               ��Ⱥ����
%����i��               ��i���˹���
%����D��               �������
%����Visual��          ��֪����
%����deta��            ӵ��������
%����trynumber��       �����̽����
%���Xinext��          ���ҵ���·��
%���flag��            ����Ƿ��ҵ����õ�·����flag=0��ʾ��Ⱥʧ�ܣ�flag=1��ʾ��Ⱥ�ɹ�
function [Xinext,flag]=AF_swarm(X,i,D,Visual,deta,trynumber)

Xi=X(i,:);                                                      %��i���˹���
N=size(X,1);                                                    %��Ⱥ��Ŀ
Yi=PathLength(D,Xi);                                            %·��Xi���ܾ���
neighbork=k_neighborhood(X,i,Visual);                           %Xi�����򼯺�
nf=size(neighbork,1);                                           %���򼯺����������
flag=0;                                                         %����Ƿ��Ⱥ�ɹ�

Xc=Center(neighbork);                                           %neighbork�����ġ�·����
if ~isempty(Xc)
    Yc=PathLength(D,Xc);                                            %·��Xc���ܾ���
    if (Yc<Yi)&&(nf/N<deta)
        Xinext=Xc;
        flag=1;
    else
        [Xinext,flag]=AF_prey(X,i,D,trynumber,Visual);
    end
else
    [Xinext,flag]=AF_prey(X,i,D,trynumber,Visual);
end
end

