%% ��ʼ����Ⱥ
%����FishNum��                         ��Ⱥ����
%����num_Citys��                       ��������
%���initFish:                         ��ʼ����Ⱥ
function initFish=AF_init(FishNum,num_Citys)
initFish=zeros(FishNum,num_Citys);
for i=1:FishNum
    initFish(i,:)=randperm(num_Citys);
end
end

