%% ����������ʹ��i���˹�����ʳ,������ʳ�ɹ�����flag Ϊ1��X��Ϊi����ʳ���״̬������flagΪ0
%����X��               ��Ⱥ����
%����i��               ��i���˹���
%����D��               �������
%����trynumber��       �����̽����
%����Visual��          ��֪����
%���Xinext��          ���ҵ���·��
%���flag��            ����Ƿ��ҵ����õ�·����flag=0��ʾ��ʳʧ�ܣ�flag=1��ʾ��ʳ�ɹ�
function [Xinext,flag]=AF_prey(X,i,D,trynumber,Visual)

Xinext=[];
Yi=PathLength(D,X(i,:));                                            %·��Xi���ܾ���
CityNum=length(X(i,:));
flag=0;                                                         %����Ƿ���ʳ������·����flag=0��ʾû��ʳ����flag=1��ʾ��ʳ�ɹ�
for j=1:trynumber
    while(1)
        DJ=floor(rand*Visual)+1;    %����ͬ���ֶ���
        if(DJ>0 && DJ<=Visual)         %����Ұ��Χ��
            break;
        end
    end
    while(1)
         S(1)=floor(rand*CityNum)+1;
         if(S(1)>1 && S(1)<=CityNum)  %�����г�����
            break;
         end
    end
    p=1;
    while(p<DJ)
       t=floor(rand*CityNum)+1;
       if(t>1&&t<=CityNum && sum(S==t)==0)
           p=p+1;
           S(p)=t;
       end
   end
   Xi=X(i,:);
   t=Xi(S(1));
   for k=1:DJ-1
       Xi(S(k))=Xi(S(k+1));
   end
   Xi(S(DJ))=t;
   YY=PathLength(D,Xi);                                            %·��Xi���ܾ���
   if YY<Yi
       Xinext=Xi;
       flag=1;
       return;
   end
end
Xinext=Xi;

end

