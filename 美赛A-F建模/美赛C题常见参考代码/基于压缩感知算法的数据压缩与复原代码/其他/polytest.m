x=input('x=');
y=input('y=');
nx=length(x);
ny=length(y);
if nx==ny
    x1=x(1);xn=x(nx);
    disp('ͨ������Ľ���ʽͼ�Σ������Ԥ������Ҫ��϶���ʽ�Ľ���...');
    disp('�밴�س�����polytool');
 pause;
 polytool(x,y,1);
 m=input('�������ʽ�Ľ���m=');
 [p,s]=polyfit(x,y,m);%�������
 disp('�������ʽ�ĸ���ϵ��');
 fprintf(1,'a=%3.16f\n',p);
 disp '�������ʽ�������Ϣs'
 disp 's'
 [y1,delta]=polyconf(p,x,s);%����ϵõ��Ķ���ʽ����������
 disp('�۲���ϵ�����');
 disp('x--------- y--------- y1------------');
 for i=1:nx
     xy=[x(i) y(i) y1(i)];
     disp(xy)
 end
 plot(x,y,'r.');
 title('����ʽ���ʵ��');
 hold on;
 xi=[x1:0.1:xn];
 yi=polyval(p,xi);
 plot(xi,yi,'k-');
 %���Ч���;��ȼ���
 disp '���������ʾ'
 Q=norm(y-y1).^2

else 
    disp '����������������������г���������ȷ������'��
 clear
end
 
 