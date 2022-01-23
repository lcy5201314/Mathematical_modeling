% ==========================���ٽ����㷨===============================
% ================���������ʹ�õ���Ҫ����==============================
% X��  ��������������� 
% l:   ���ֵ��Ӽ���Ŀ
% L:   ˮƽ��Ŀ
% Xp:  �ڵ�p��Ӧ�������Ӽ�
% Mp:   ����ľ�ֵ
% Rp:   ��Mp��Xi����Զ����
% =============������������ʹ�õ���Ҫ����=================================
% CurL:  ��ǰˮƽ
% p:     ��ǰ���
% CurTable:  ��ǰĿ¼���е���������
% CurPinT:   �ڵ�ǰĿ¼���е����������
% RpCur:     ��ǰĿ¼���н��p��Ӧ��Rp
% x:        ��������
% =====================================================================
% ʵ�������������㷨�ھ������֮�󣬽������������ٶȵ�ȷ��һ��Ľ��ڷ�����
% ��ʱ���ھ���Ҫ���Ĵ�����ʱ�䣬������ٶȲ���һ��Ľ��ڷ���
% =====================================================================
%   Copyright Wang Chuanting.
%   $Revision: 1.0 $  $Date: 2008/05/09 09:40:34 $
% ======================================================================

% ======================================================================
% ���������Ƚ��о��࣭����
clear,close all;
% tic
X = [randn(200,2)+ones(200,2);...
     randn(200,2)-2*ones(200,2);...
     randn(200,2)+4*ones(200,2);];
% ������ÿ��ˮƽ������Ϊl���Ӽ�������
[row,col]=size(X);
all_idx=0;L=3;l=3;    
%�����ܽڵ����Ŀ
for i=1:L
    all_idx=all_idx+l^i;
end
Xp=cell(all_idx,1);
Mp=zeros(all_idx,col);
Rp=zeros(all_idx,1);
p=1;
for i=1:L
   if i==1
       [IDX,C,sumd,D] = kmeans(X,l);
        for j=1:l
            Xp(p)={X((IDX==j),:)};
            Mp(p,:)=C(j,:);
            Rp(p)=max(D((IDX==j),j));
            p=p+1;
        end   
   else
       endk=p-1;begink=endk-l^(i-1)+1;       
       for k=begink:endk           
           [IDX,C,sumd,D] = kmeans(Xp{k,1},l);
           X1=Xp{k,1};
           for j=1:l                
                Xp(p)={X1((IDX==j),:)};
                Mp(p,:)=C(j,:);
                Rp(p)=max(D((IDX==j),j));
                p=p+1;
           end   
       end
   end
end
% ====================================================================
% ����������������������
tic
x=randn(1,2);%��������
B=inf;CurL=1;p=0;TT=1;
while TT==1 %����2
    Xcurp=cell(1);
    CurTable=cell(l,1);
    CurPinT=zeros(l,1);
    Dx=zeros(l,1);
    RpCur=zeros(l,1);
    %��ǰ�ڵ��ֱ�Ӻ�̷���Ŀ¼��   
    for i=1:l   
        CurTable(i,1)=Xp(i+p*l,1);
        CurPinT(i)=i+p*l;
        Dx(i)=norm(x-Mp(i+p*l,:))^2;
        RpCur(i)=Rp(i+p*l);
    end    
    
    while 1 %����3
        [rowT,colT]=size(CurTable);
        for i=1:rowT                   
            if Dx(i)>B+RpCur(i)+eps%��Ŀ¼����ȥ����ǰ�ڵ�p
                CurTable(i,:)=[];
                CurPinT(i)=[];
                Dx(i)=[];
                RpCur(i)=[];
                break;
            end
        end
        [CurRowT,CurColT]=size(CurTable);
        if CurRowT==0
           CurL=CurL-1;p=floor((p-1)/3);
           if CurL==0
              TT=0; break;  
           else
               %ת����3
           end
        elseif CurRowT>0
            [Dxx,Dxind]=sort(Dx,'ascend');
            p1=CurPinT(Dxind(1));
            p=p1;
            %�ӵ�ǰĿ¼��ȥ��p1
            for j=1:CurRowT
                if CurPinT(j)==p1
                    Xcurp(1,1)=CurTable(j,1);
                    CurTable(j,:)=[];
                    CurPinT(j)=[];
                    CurD=Dx(j);%��¼D(x,Mp)
                    Dx(j)=[];
                    RpCur(j)=[];                    
                    break;
                end
            end
            if CurL==L
                XcurpMat=cell2mat(Xcurp);
                [CurpRow,CurpCol]=size(XcurpMat);
                CurpMean=Mp(p,:);
                for k=1:CurpRow
                    Dxi=norm((XcurpMat(k,:)-CurpMean))^2;
                    if CurD>Dxi+B+eps
                        
                    else
                        Dxxi=norm((x-XcurpMat(k,:)))^2;
                        if Dxxi<B+eps
                            B=Dxxi;Xnn=XcurpMat(k,:);
                        end
                    end
                end
            else
                CurL=CurL+1;
                break;
            end
        end
    end
end
B,Xnn,NN=find(X(:,1)==Xnn(1))
time1=toc
% ====================================================================
figure, plot(X(1:200,1),X(1:200,2),'m.')
hold on,plot(X(201:400,1),X(201:400,2),'b.')
hold on,plot(X(401:600,1),X(401:600,2),'g.')   
hold on,plot(Xnn(1),Xnn(2),'kx ','MarkerSize',10,'LineWidth',2)
hold on,plot(x(1),x(2),'r+','MarkerSize',10,'LineWidth',2)
legend('Cluster 1','Cluster 2','Cluster 3','NN','x','Location','NW')


