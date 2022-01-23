x=input('x=');
y=input('y=');
nx=length(x);
ny=length(y);
if nx==ny
    x1=x(1);xn=x(nx);
    disp('通过下面的交互式图形，你可以预估计下要拟合多项式的阶数...');
    disp('请按回车键打开polytool');
 pause;
 polytool(x,y,1);
 m=input('输入多项式的阶数m=');
 [p,s]=polyfit(x,y,m);%曲线拟合
 disp('输出多项式的各项系数');
 fprintf(1,'a=%3.16f\n',p);
 disp '输出多项式的相关信息s'
 disp 's'
 [y1,delta]=polyconf(p,x,s);%用拟合得到的多项式计算纵坐标
 disp('观测拟合的数据');
 disp('x--------- y--------- y1------------');
 for i=1:nx
     xy=[x(i) y(i) y1(i)];
     disp(xy)
 end
 plot(x,y,'r.');
 title('多项式拟合实验');
 hold on;
 xi=[x1:0.1:xn];
 yi=polyval(p,xi);
 plot(xi,yi,'k-');
 %拟合效果和精度检验
 disp '均方误差显示'
 Q=norm(y-y1).^2

else 
    disp '输入的数据有误，请重新运行程序并输入正确的数据'；
 clear
end
 
 