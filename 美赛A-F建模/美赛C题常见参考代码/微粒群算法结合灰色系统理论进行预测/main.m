%���ñ��ż��㽭��ũ������Ԥ��
%2010-4-4
clear all
clc
%N=3
k=0;
n=2;
m=100;
c1=2;
c2=2;
c=1;
w=0.729;
iter=1;
itermax=1500;%gai1
wmax=1.4;
wmin=0.35;
vmax=2;

x=-1+2*rand(m,n); %ԭ���ǲ���25��2�е�-1��1��Χ�������
v=2*rand(m,n);%����25��2�е�0��2��Χ�������

fid = fopen('xb.txt', 'r');
data = fscanf(fid, '%f',12);   % ��1.txt�е�21��������˫���ȶ���
fclose(fid);
 w=wmax-iter*(wmax-wmin)/itermax; 
for i=1:m
    [f(i) t]=huise(data,x(i,:));%f(i)ΪԤ��ֵ��ԭֵ֮�������ۼ�ֵ��tΪԤ��ֵ
end

pbx=x;
pbf=f;
[gbf i]=min(pbf);
gbx=pbx(i,:);
for i=1:m
    v(i,:)=w*v(i,:)+c1*rand*(pbx(i,:)-x(i,:))+c2*rand*(gbx-x(i,:));%�����ٶ�
    for j=1:n 
        if v(i,j)>vmax
            v(i,j)=vmax;
        elseif v(i,j)<-vmax
            v(i,j)=-vmax; 
        end 
    end  %��ÿһ�е��ٶ�����������vmax֮��
    x(i,:)=x(i,:)+v(i,:); %����λ��
end   %Ӧ���ǳ�ʼ����


while  abs(gbf)>0.0001 && iter<=1500%gai2
     
    for i=1:m
        [f(i) t]=huise(data,x(i,:));
    end
    for i=1:m
        if f(i)<pbf(i)
            pbf(i)=f(i);
            pbx(i,:)=x(i,:);
        end
    end       %�ھ��ȼ�����δ��Ҫ�����ٴ�Ԥ��
    [gbf i]=min(pbf);
    gbx=pbx(i,:);
    
 %   if bianbie(f)<c
 %    [gbx gbf]=hundun(data,gbx);
 %   end
 
w=wmax-iter*(wmax-wmin)/itermax;

    for i=1:m
         v(i,:)=w*v(i,:)+c1*rand*(pbx(i,:)-x(i,:))+c2*rand*(gbx-x(i,:));
         for j=1:n
             if v(i,j)>vmax
            v(i,j)=vmax;
        elseif v(i,j)<-vmax
            v(i,j)=-vmax;
        end
         end
     x(i,:)=x(i,:)+v(i,:);
    end
iter=iter+1;
end  %����while����

a=gbx(1) %
b=gbx(2) %
f=abs(gbf)%���������ֵ֮��
yuan=data
yuce=t'
T=length(data);
pw=0;
for i=1:T
xdwc(1:T,1)=((yuan(1:T,1)-yuce(1:T,1))./yuan)*100 ;%������
pw=pw+xdwc(i);
end
xdwc
pw=pw/T    %ƽ��������

figure(1)
plotljz(data,t(1:T));

% %t;

    