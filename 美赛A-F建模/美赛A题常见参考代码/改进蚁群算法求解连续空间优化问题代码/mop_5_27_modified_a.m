%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Ant Colony System Code                          %
%                      Date 5/27/2006                              %
%                       Theobald.Zou                               %
%                      Original  Code                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
tic;
%��ʼ��������
Ant=50;%���Ϲ�ģ
ECHO=50;%��������
%F=(X1.-1)^2+(X2.-2.2)^2+1;%���Ż�����
step=0.1*rand(1);%�ֲ�����ʱ�Ĳ���
temp=[0,0];
%�������䳤��
start1=0;
end1=2;
start2=1;
end2=3;
Len1=(end1-start1)/Ant;
Len2=(end2-start2)/Ant;
%P = 0.2;
%��ʼ������λ��
subplot(2,2,1);
for i=1:Ant
X(i,1)=(start1+(end1-start1)*rand(1));
X(i,2)=(start2+(end2-start2)*rand(1));
T0(i)=exp(-((X(i,1)-1)^2+(X(i,2)-2.2)^2+1));%��ʼ��Ϣ��,�溯��ֵ��,��Ϣ��Ũ��С,��֮��Ȼ
plot(X(i,1),X(i,2),'k.')
hold on;
title('(a)');
xlabel('X1');
ylabel('X2');
grid on;
end;
%���˳�ʼ�����
for Echo=1:ECHO    %��ʼѰ��    
%P0��������,P0Ϊȫ��ת��ѡ������
a1=0.9;
b1=(1/ECHO)*2*log(1/2);
f1=a1*exp(b1*Echo);
a2=0.225;
b2=(1/ECHO)*2*log(2);
f2=a2*exp(b2*Echo);
if  Echo<=(ECHO/2)
    P0=f1;
else
    P0=f2;
end;
%P�������壬PΪ��Ϣ������ϵ��
a3=0.1;
b3=(1/ECHO).*log(9);
P=a3*exp(b3*Echo);
lamda=0.10+(0.14-0.1)*rand(1);%ȫ��ת�Ʋ�������
Wmax=1.0+(1.4-1.0)*rand(1);%�������²�������
Wmin=0.2+(0.8-0.2)*rand(1);%�������²�������
%Ѱ�ҳ�ʼ����ֵ
   T_Best=T0(1);
    for j=1:Ant
       if  T0(j)>=T_Best
           T_Best=T0(j);
           BestIndex=j;
       end;    
    end;   
    W=Wmax-(Wmax-Wmin)*(Echo/ECHO);  %�ֲ������������²���
    for j_g=1:Ant  %ȫ��ת�Ƹ�����ȡ,������������λ�ò���bestindexʱ
        if j_g~=BestIndex
            r=T0(BestIndex)-T0(j_g);
            Prob(j_g)=exp(r)/exp(T0(BestIndex));
        else   %��j_g=BestIndex��ʱ����оֲ�����
            if rand(1)<0.5
                temp(1,1)=X(BestIndex,1)+W*step;
                temp(1,2)=X(BestIndex,2)+W*step;
            else
                temp(1,1)=X(BestIndex,1)-W*step;
                temp(1,2)=X(BestIndex,2)-W*step;
            end;
            Prob(j_g)=0;%bestindex�����ϲ�����ȫ��ת��
        end;
            X1_T=temp(1,1);
            X2_T=temp(1,2);
            X1_B=X(BestIndex,1);
            X2_B=X(BestIndex,2);
            F1_T=(X1_T-1).^2+(X2_T-2.2).^2+1;
            F1_B=(X1_B-1).^2+(X2_B-2.2).^2+1;
            if exp(-F1_T)>exp(-F1_B)
               X(BestIndex,1)=temp(1,1);
               X(BestIndex,2)=temp(1,2);
            end;
    end;  
    for j_g_tr=1:Ant
        if Prob(j_g_tr)<P0
            X(j_g_tr,1)=X(j_g_tr,1)+lamda*(X(BestIndex,1)-X(j_g_tr,1));%Xi=Xi+lamda*(Xbest-Xi)
            X(j_g_tr,2)=X(j_g_tr,2)+lamda*(X(BestIndex,2)-X(j_g_tr,2));%Xi=Xi+lamda*(Xbest-Xi)
        else
            X(j_g_tr,1)=X(j_g_tr,1)+((-1)+2*rand(1))*Len1;%Xi=Xi+rand(-1,1)*Len1
            X(j_g_tr,2)=X(j_g_tr,2)+((-1)+2*rand(1))*Len2;%Xi=Xi+rand(-1,1)*Len2
        end;
    end;
%��Ϣ�ظ���
    for t_t=1:Ant
         T0(t_t)=(1-P)*T0(t_t)+(exp(-(X(t_t,1)-1).^2+(X(t_t,2)-2.2).^2+1));
    end;
    if Echo==round(ECHO/3)%��������1/3ʱ����ɫ���ʾ���ϵķֲ�λ��
        subplot(2,2,2);
        for i_draw1=1:Ant
        plot(X(i_draw1,1),X(i_draw1,2),'g.')
        axis([0 2 1 3]);
        hold on;
        title('(b)');
        xlabel('X1');
        ylabel('X2');
        end;
        grid on;
    end;
    [c_iter,i_iter]=max(T0);  %��ȡÿ��ȫ�����Ž�
    minpoint_iter=[X(i_iter,1),X(i_iter,2)];
    minvalue_iter=(X(i_iter,1)-1).^2+(X(i_iter,2)-2.2).^2+1;
    min_local(Echo)=minvalue_iter;%����ÿ���ֲ����Ž�
   %��ÿ��ȫ�����Ž�浽min_global������
   if Echo >= 2
    if min_local(Echo)<min_global(Echo-1)
       min_global(Echo)=min_local(Echo);
    else
        min_global(Echo)=min_global(Echo-1);
    end;
   else
    min_global(Echo)=minvalue_iter;
   end;
 end;%ECHOѭ������
  subplot(2,2,3);
  for i_draw3=1:Ant
  plot(X(i_draw3,1),X(i_draw3,2),'r.')%���������ú�ɫ���ʾ���ϵķֲ�λ��
  axis([0 2 1 3]);
  hold on;
  title('(c)');
  xlabel('X1');
  ylabel('X2');
  end;
  grid on;
  subplot(2,2,4);
  min_global=min_global';
  index(:,1)=1:ECHO;
  plot(index(:,1), min_global(:,1),'b-')
  hold on;
  title('(d)');
  xlabel('iteration');
  ylabel('f(x)');
  grid on;
  [c_max,i_max]=max(T0);  
  minpoint=[X(i_max,1),X(i_max,2)]
  minvalue=(X(i_max,1)-1).^2+(X(i_max,2)-2.2).^2+1
  runtime=toc