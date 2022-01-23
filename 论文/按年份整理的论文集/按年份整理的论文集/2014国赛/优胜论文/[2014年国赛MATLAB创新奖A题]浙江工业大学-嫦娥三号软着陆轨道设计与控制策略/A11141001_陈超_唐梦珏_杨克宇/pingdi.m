% �ҳ���ͼ��ƽ��
clc,clear all
f1=imread('����3 ��2400m�������ָ߳�ͼ.tif');
f2=imread('����4 ������100m�������ָ߳�ͼ.tif');
[n,m]=size(f1);
count=zeros(max(max(f1))+1,1);
for i=1:n
    for j=1:m
        count(f1(i,j)+1) = count(f1(i,j)+1)+1;
    end
end
[m,h]=max(count);
[n,m] = size(f1);
for i=1:n
    for j=1:m
        if f1(i,j)>h-1+4 || f1(i,j)<h-1-4
            f1(i,j)=-1;
        end
    end
end
x=1:m;
y=1:n;
mesh(x', y', double(f1));

[n,m]=size(f2);
count=zeros(max(max(f2))+1,1);
for i=1:n
    for j=1:m
        count(f2(i,j)+1) = count(f2(i,j)+1)+1;
    end
end
[m,h]=max(count);
[n,m] = size(f2);
for i=1:n
    for j=1:m
        if f2(i,j)>h-1+4 || f2(i,j)<h-1-4
            f2(i,j)=-1;
        end
    end
end
x=1:m;
y=1:n;
figure,mesh(x', y', double(f2));
save('pingdif4.mat','f1','f2');