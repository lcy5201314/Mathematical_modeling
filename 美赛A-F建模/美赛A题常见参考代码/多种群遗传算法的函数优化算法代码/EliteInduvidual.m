%% �˹�ѡ������
%����Chrom��           cell���ͣ�ÿ��Ԫ����ԪΪһ����Ⱥ�ı��루����ǰ��
%����ObjV��            cell���ͣ�ÿ��Ԫ��Ϊһ����Ⱥ���и����Ŀ��ֵ������ǰ��
%����MaxObjV��         Double���ͣ�ÿ����Ⱥ��ǰ���Ÿ����Ŀ��ֵ��ѡ��ǰ��
%����MaxChrom��        Double���ͣ�ÿ����Ⱥ��ǰ���Ÿ���ı��루ѡ��ǰ��
%����MaxObjV��         Double���ͣ�ÿ����Ⱥ��ǰ���Ÿ����Ŀ��ֵ��ѡ���
%����MaxChrom��        Double���ͣ�ÿ����Ⱥ��ǰ���Ÿ���ı��루ѡ���

function [MaxObjV,MaxChrom]=EliteInduvidual(Chrom,ObjV,MaxObjV,MaxChrom)

MP=length(Chrom);  %��Ⱥ��
for i=1:MP
    [MaxO,maxI]=max(ObjV{i});   %�ҳ���i��Ⱥ�����Ÿ���
    if MaxO>MaxObjV(i)
        MaxObjV(i)=MaxO;         %��¼����Ⱥ�ľ�������
        MaxChrom(i,:)=Chrom{i}(maxI,:);  %��¼����Ⱥ��������ı���
    end
end
