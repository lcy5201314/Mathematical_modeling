function [Chrom,ObjV]=immigrant(Chrom,ObjV)
%% ��������
MP=length(Chrom);
for i=1:MP
    [MaxO,maxI]=max(ObjV{i});  % �ҳ���i��Ⱥ�����ŵĸ���
    next_i=i+1;                % Ŀ����Ⱥ����������У�
    if next_i>MP;next_i=mod(next_i,MP);end
    [MinO,minI]=min(ObjV{next_i});          %  �ҳ�Ŀ����Ⱥ�����ӵĸ���
    %% Ŀ����Ⱥ���Ӹ����滻ΪԴ��Ⱥ���Ÿ���
    Chrom{next_i}(minI,:)=Chrom{i}(maxI,:);
    ObjV{next_i}(minI)=ObjV{i}(maxI);
end

