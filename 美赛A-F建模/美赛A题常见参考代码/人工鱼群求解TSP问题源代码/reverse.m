%% reverse���� ��������sequence������i��kλ��֮�������������
% ���룺Xi                         ��ʼ��������
% ���룺D                          �������
% �����R                          ������������������
function [R,flagi]= reverse(Xi,D)

num_Citys=length(Xi);                                           %������Ŀ
Yi=PathLength(D,Xi);                                            %·��Xi���ܾ���
flagi=0;
flagk=0;
for i=1:num_Citys-1
    for k=i+1:num_Citys
        XRev=Xi;
        XRev(i:k)=Xi(k:-1:i);
        YRev=PathLength(D,XRev);                                   %·��XRev���ܾ���
        if YRev<Yi
            R=XRev;
            flagk=1;
            break
        end
    end
    if flagk==1
        flagi=1;
        break
    end
end

if flagi==0
    R=randperm(num_Citys);
end

end

