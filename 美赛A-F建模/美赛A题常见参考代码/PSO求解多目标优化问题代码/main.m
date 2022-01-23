%% �ú�����ʾ��Ŀ��perota�Ż�����
%��ջ���
clc
clear

load data


%% ��ʼ����
objnum=size(P,1); %������Ʒ����
weight=92;        %����������

%��ʼ������
Dim=5;     %����ά��
xSize=50;  %��Ⱥ����
MaxIt=200; %��������
c1=0.8;    %�㷨����
c2=0.8;    %�㷨���� 
wmax=1.2;  %��������
wmin=0.1;  %��������

x=unidrnd(4,xSize,Dim);  %���ӳ�ʼ��
v=zeros(xSize,Dim);      %�ٶȳ�ʼ��

xbest=x;           %�������ֵ
gbest=x(1,:);      %����Ⱥ���λ��

% ������Ӧ��ֵ 
px=zeros(1,xSize);   %���Ӽ�ֵĿ��
rx=zeros(1,xSize);   %�������Ŀ��
cx=zeros(1,xSize);   %����Լ��

% ����ֵ��ʼ��
pxbest=zeros(1,xSize); %�������ż�ֵĿ��
rxbest=zeros(1,xSize); %�����������Ŀ��
cxbest=zeros(1,xSize);  %��¼����������Լ��

% ��һ�ε�ֵ
pxPrior=zeros(1,xSize);%���Ӽ�ֵĿ��
rxPrior=zeros(1,xSize);%�������Ŀ��
cxPrior=zeros(1,xSize);%��¼����������Լ��

%�����ʼĿ������
for i=1:xSize
    for j=1:Dim %�������
        px(i) = px(i)+P(x(i,j),j);  %���Ӽ�ֵ
        rx(i) = rx(i)+R(x(i,j),j);  %�������
        cx(i) = cx(i)+C(x(i,j),j);  %��������
    end
end
% ��������λ��
pxbest=px;
rxbest=rx;
cxbest=cx;

%% ��ʼɸѡ���ӽ�
flj=[];
fljx=[];
fljNum=0;
%����ʵ����Ⱦ���
tol=1e-7;
for i=1:xSize
    flag=0;  %֧���־
    for j=1:xSize  
        if j~=i
            if ((px(i)<px(j)) &&  (rx(i)>rx(j))) ||((abs(px(i)-px(j))<tol)...
                    &&  (rx(i)>rx(j)))||((px(i)<px(j)) &&  (abs(rx(i)-rx(j))<tol)) || (cx(i)>weight) 
                flag=1;
                break;
            end
        end
    end
    
    %�ж����ޱ�֧��
    if flag==0
        fljNum=fljNum+1;
        % ��¼���ӽ�
        flj(fljNum,1)=px(i);
        flj(fljNum,2)=rx(i);
        flj(fljNum,3)=cx(i);
        % ���ӽ�λ��
        fljx(fljNum,:)=x(i,:); 
    end
end

%% ѭ������
for iter=1:MaxIt
    
    % Ȩֵ����
    w=wmax-(wmax-wmin)*iter/MaxIt;
     
    %�ӷ��ӽ���ѡ��������Ϊȫ�����Ž�
    s=size(fljx,1);       
    index=randi(s,1,1);  
    gbest=fljx(index,:);

    %% Ⱥ�����
    for i=1:xSize
        %�ٶȸ���
        v(i,:)=w*v(i,:)+c1*rand(1,1)*(xbest(i,:)-x(i,:))+c2*rand(1,1)*(gbest-x(i,:));
        
        %λ�ø���
        x(i,:)=x(i,:)+v(i,:);
        x(i,:) = rem(x(i,:),objnum)/double(objnum);
        index1=find(x(i,:)<=0);
        if ~isempty(index1)
            x(i,index1)=rand(size(index1));
        end
        x(i,:)=ceil(4*x(i,:));        
    end
    
    %% ���������Ӧ��
    pxPrior(:)=0;
    rxPrior(:)=0;
    cxPrior(:)=0;
    for i=1:xSize
        for j=1:Dim %�������
            pxPrior(i) = pxPrior(i)+P(x(i,j),j);  %��������i ��ֵ
            rxPrior(i) = rxPrior(i)+R(x(i,j),j);  %��������i ���
            cxPrior(i) = cxPrior(i)+C(x(i,j),j);  %��������i ����
        end
    end
    
    %% ����������ʷ���
    for i=1:xSize
        %���ڵ�֧��ԭ�еģ����ԭ�е�
         if ((px(i)<pxPrior(i)) &&  (rx(i)>rxPrior(i))) ||((abs(px(i)-pxPrior(i))<tol)...
                 &&  (rx(i)>rxPrior(i)))||((px(i)<pxPrior(i)) &&  (abs(rx(i)-rxPrior(i))<tol)) || (cx(i)>weight) 
                xbest(i,:)=x(i,:);%û�м�¼Ŀ��ֵ
                pxbest(i)=pxPrior(i);
                rxbest(i)=rxPrior(i);
                cxbest(i)=cxPrior(i);
          end
        
        %�˴˲���֧�䣬�������
        if ~( ((px(i)<pxPrior(i)) &&  (rx(i)>rxPrior(i))) ||((abs(px(i)-pxPrior(i))<tol)...
                &&  (rx(i)>rxPrior(i)))||((px(i)<pxPrior(i)) &&  (abs(rx(i)-rxPrior(i))<tol)) || (cx(i)>weight) )...
                &&  ~( ((pxPrior(i)<px(i)) &&  (rxPrior(i)>rx(i))) ||((abs(pxPrior(i)-px(i))<tol) &&  (rxPrior(i)>rx(i)))...
                ||((pxPrior(i)<px(i)) &&  (abs(rxPrior(i)-rx(i))<tol)) || (cxPrior(i)>weight) )
            if rand(1,1)<0.5
                xbest(i,:)=x(i,:);
                  pxbest(i)=pxPrior(i);
                  rxbest(i)=rxPrior(i);
                  cxbest(i)=cxPrior(i);
            end
        end
    end

    %% ���·��ӽ⼯��
    px=pxPrior;
    rx=rxPrior;
    cx=cxPrior;
    %�����������ӽ⼯��
    s=size(flj,1);%Ŀǰ���ӽ⼯����Ԫ�ظ���
   
    %�Ƚ����ӽ⼯�Ϻ�xbest�ϲ�
    pppx=zeros(1,s+xSize);
    rrrx=zeros(1,s+xSize);
    cccx=zeros(1,s+xSize);
    pppx(1:xSize)=pxbest;
    pppx(xSize+1:end)=flj(:,1)';
    rrrx(1:xSize)=rxbest;
    rrrx(xSize+1:end)=flj(:,2)';
    cccx(1:xSize)=cxbest;
    cccx(xSize+1:end)=flj(:,3)';
    xxbest=zeros(s+xSize,Dim);
    xxbest(1:xSize,:)=xbest;
    xxbest(xSize+1:end,:)=fljx;
   
    %ɸѡ���ӽ�
    flj=[];
    fljx=[];
    k=0;
    tol=1e-7;
    for i=1:xSize+s
        flag=0;%û�б�֧��
        %�жϸõ��Ƿ����
        for j=1:xSize+s 
            if j~=i
                if ((pppx(i)<pppx(j)) &&  (rrrx(i)>rrrx(j))) ||((abs(pppx(i)-pppx(j))<tol) ...
                        &&  (rrrx(i)>rrrx(j)))||((pppx(i)<pppx(j)) &&  (abs(rrrx(i)-rrrx(j))<tol)) ...
                        || (cccx(i)>weight) %��һ�α�֧��
                    flag=1;
                    break;
                end
            end
        end

        %�ж����ޱ�֧��
        if flag==0
            k=k+1;
            flj(k,1)=pppx(i);
            flj(k,2)=rrrx(i);
            flj(k,3)=cccx(i);%��¼���ӽ�
            fljx(k,:)=xxbest(i,:);%���ӽ�λ��
        end
    end
    
    %ȥ���ظ�����
    repflag=0;   %�ظ���־
    k=1;         %��ͬ���ӽ�������
    flj2=[];     %�洢��ͬ���ӽ�
    fljx2=[];    %�洢��ͬ���ӽ�����λ��
    flj2(k,:)=flj(1,:);
    fljx2(k,:)=fljx(1,:);
    for j=2:size(flj,1)
        repflag=0;  %�ظ���־
        for i=1:size(flj2,1)
            result=(fljx(j,:)==fljx2(i,:));
            if length(find(result==1))==Dim
                repflag=1;%���ظ�
            end
        end
        %���Ӳ�ͬ���洢
        if repflag==0 
            k=k+1;
            flj2(k,:)=flj(j,:);
            fljx2(k,:)=fljx(j,:);
        end
        
    end
    
    %���ӽ����
    flj=flj2;
    fljx=fljx2;

end

%���Ʒ��ӽ�ֲ�
plot(flj(:,1),flj(:,2),'o') 
xlabel('P')
ylabel('R')
title('���շ��ӽ���Ŀ��ռ�ֲ�')
disp('���ӽ�flj����������ΪP��R��C')