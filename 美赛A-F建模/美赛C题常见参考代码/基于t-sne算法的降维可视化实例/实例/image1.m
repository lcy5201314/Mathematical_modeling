clc
clear
w='D:\image\2\2_-1.bmp'%读取文件地址
k=1
for i=0:9
    for j=1:500
        a=num2str(i);
        b=num2str(j);
        w1=strrep(w,'2',a);
        w2=strrep(w1,'-1',b);
        x=imread(w2);
        y(k,:)=x(:)';
        z(k)=i;
        k=k+1;
    end
end
y1=double(y);


 x=tsne1(y1);
 gscatter(x(:,1),x(:,2),z');

 x=tsne2(y1,'Algorithm','exact' );
figure;
 gscatter(x(:,1),x(:,2),z');
