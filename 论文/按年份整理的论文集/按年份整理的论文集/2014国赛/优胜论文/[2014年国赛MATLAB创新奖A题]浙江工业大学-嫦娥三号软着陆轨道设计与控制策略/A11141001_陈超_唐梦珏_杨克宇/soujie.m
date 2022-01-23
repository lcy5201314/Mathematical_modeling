clc,clear all;
% tic
global F;                   % ����
global beta;                % ���ĽǶ�
global mu;                  % ������������
global c;                   % �ȳ�
mu=4.903737416799999e+12;   % ������������
c=2940;                     % �ȳ�
vf=zeros(401,360); 
for li=1500:7500
    F=li;
    for ja=1:3600
        beta=pi/180*ja*0.1;             % ���ĽǶ�
        [t,x]=ode45(@func,[0:1:400],[1753000;0;0;1692/1753000;2400]);
        h=x(:,1)-1738000;
        [n,m]=size(h);
        for i=1:n
            vf(i,ja)=sqrt((x(i,4)*(h(i)+1737000))^2+x(i,2)^2);
            if h(i)>2900 & h(i)<3100
                if vf(i,ja)>47 & vf(i,ja)<67
                    ja
                end
            end
        end
    end
end



% toc
%
% figure,plot(t,h);title('h');
% figure,plot(t,x(:,2));title('v');
% figure,plot(t,x(:,3));title('cta');
% figure,plot(t,x(:,4));title('omg');
% figure,plot(t,x(:,5));title('m');