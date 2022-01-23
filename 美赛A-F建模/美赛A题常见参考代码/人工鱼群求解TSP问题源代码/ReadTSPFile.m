function [ n_citys,city_position ] = ReadTSPFile( filename )
%READTSPFILE ��ȡTSP�ļ���Ϣ
%   filename :TSP�ļ���
%   n_city  : ���и���
%   city_position ��������
 fid = fopen(filename,'rt');  %���ı�ֻ����ʽ���ļ�
 if(fid<=0)
     disp('�ļ���ʧ�ܣ�')
     return;
 end
 location=[];A=[1 2];
 tline = fgetl(fid);%��ȡ�ļ���һ��
while ischar(tline)
    if(strcmp(tline,'NODE_COORD_SECTION'))
        while ~isempty(A)
            A=fscanf(fid,'%f',[3,1]);%��ȡ�ڵ��������ݣ�ÿ�ζ�ȡһ��֮���ļ�ָ����Զ�ָ����һ��
            if isempty(A)
                break;
            end
            location=[location;A(2:3)'];%���ڵ�����浽location��
        end
    end
    tline = fgetl(fid);
    if strcmp(tline,'EOF')   %�ж��ļ��Ƿ����
        break;
    end
end
[m,n]=size(location);
n_citys=m;
city_position =location;
 fclose(fid);
end
