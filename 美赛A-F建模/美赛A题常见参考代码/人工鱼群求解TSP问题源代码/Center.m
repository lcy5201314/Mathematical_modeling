%% ���������·���ġ����ġ�·��
%����neighbork��           ����
%���center_route��        ����·��
function center_route=Center(neighbork)
num_Citys=size(neighbork,2);        %������Ŀ
XC=[];
for j=1:num_Citys
    tJ=neighbork(:,j);%���γ����ҳ��������ĳ��У��Ӷ�ȷ����������
    while(~isempty(tJ))
        fre=[];
        for k=1:length(tJ)
            fre(k)=sum(tJ==tJ(k));
        end
        [p q]=max(fre);
        if(j==1)
            break;%����whileѭ��
        elseif(sum(XC==tJ(q))==0)
            break;
        end
        tJ(tJ==tJ(q))=[];%���γ��򲢷Ƕ��һ�٣�Ϊ���Ƿ�ֹ�ظ���ͬһ������
    end
    if(isempty(tJ))
        while(1)
            b=floor(rand*(num_Citys+1));
            if(b>0 && b<=num_Citys && sum(XC==b)==0)
                XC(j)=b;
                break;
            end
        end
    else
        XC(j)=tJ(q);
    end
end
center_route=XC;
end

