%% ��������·��֮��ġ����롱
%���룺 route1 route2
%����� ������·��֮��ġ����롱
%���ӣ�route1=[3 5 4 1 2 6]  
%      route2=[1 2 4 3 6 5]
%������·��֮��ġ����롱����Ϊ��Ӧλ�����������б�Ų�ͬ����Ŀ֮��
%route1��route2ֻ�е�3��Ԫ����ͬ��������4�����ԡ����롱Ϊ5
function Dis=distance(route1,route2)
n=length(route1);           %������Ŀ
Dis=0;
for i=1:n
    %ֻ�ж�Ӧλ����Ԫ�ز�ͬ�������롱�Ż��1
    if route1(i)~=route2(i)
        Dis=Dis+1;
    end
end

end

