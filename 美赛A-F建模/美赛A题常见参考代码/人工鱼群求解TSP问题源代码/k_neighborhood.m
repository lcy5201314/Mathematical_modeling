%% �ҳ�һ��·����k����
%����i��                           	��ѡ����ı��
%����Visual��                          ���·�������롱
%����X��                               ��Ⱥ����
%���neighbork��                       ��route�����롱Ϊk�����򼯺�
function neighbork=k_neighborhood(X,i,Visual)
neighbork=[];                           %��ʼ����Ϊ��
FishNum=size(X,1);                      %��Ⱥ��Ŀ
route=X(i,:);                           
for j=1:FishNum        
    R=X(j,:);
    Dis=distance(route,R);              %��������·��֮��ġ����롱
    if (Dis<Visual)&&(i~=j)
        neighbork(end+1,:)=R;           %ֻ�С����롱С��k�����ܳ�Ϊ�����е��� 
    end 
end

end

