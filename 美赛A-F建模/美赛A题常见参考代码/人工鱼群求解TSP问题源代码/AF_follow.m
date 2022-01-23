%% ׷β��Ϊ
%����X��               ��Ⱥ����
%����i��               ��i���˹���
%����D��               �������
%����Visual��          ��֪����
%����deta��            ӵ��������
%����trynumber��       �����̽����
%���Xinext��          ���ҵ���·��
%���flag��            ����Ƿ��ҵ����õ�·����flag=0��ʾ׷βʧ�ܣ�flag=1��ʾ׷β�ɹ�
function [Xinext,flag]=AF_follow(X,i,D,Visual,deta,trynumber)
Xi=X(i,:);                                                      %��i���˹���
N=size(X,1);                                                    %��Ⱥ��Ŀ
Yi=PathLength(D,Xi);                                            %·��Xi���ܾ���
neighbork=k_neighborhood(X,i,Visual);                           %Xi�����򼯺�
nf=size(neighbork,1);                                           %���򼯺����������
flag=0;                                                         %����Ƿ�׷β�ɹ�������·���Ƿ��ԭ��·���ܾ������

Y=AF_foodconsistence(neighbork,D);                              %�����и���·�����ܾ���
[Ymin,minIndex]=min(Y);                                         %�ҳ����򼯺����ܾ�����С������·��
Xmin=neighbork(minIndex,:);                                     %���򼯺����ܾ�����С������·��

if ~isempty(Ymin)
    if (Ymin<Yi) && (nf/N<deta)
        Xinext=Xmin;
        flag=1;
    else
        [Xinext,flag]=AF_prey(X,i,D,trynumber,Visual);
    end
else
    [Xinext,flag]=AF_prey(X,i,D,trynumber,Visual);
end

end

