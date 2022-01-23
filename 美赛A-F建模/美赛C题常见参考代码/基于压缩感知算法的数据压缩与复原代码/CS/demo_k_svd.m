
%%
%k-svd分解
%an algorithem for designing of overcomplete dictionaries for sparse  representations
%preparetion
clc;
clear;
sig_len=20;
atm_num=50;
dic=randn(sig_len,atm_num);
dic=dic./repmat(sqrt(sum(dic.^2,1)),[sig_len,1]);
[M,N]=size(dic)
Y_std=[];
for i=1:3:30
    for j=2:3:30
        for k=3:3:30
            y=2*dic(:,i)+3*dic(:,j)+4*dic(:,k);
            Y_std=[Y_std,y];
        end
    end
end
Dic_est=Y_std(:,1:atm_num);

Dic_est=Dic_est./repmat(sqrt(sum(Dic_est.^2,1)),[sig_len,1]);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%k_svd分解

Y=Y_std(:,atm_num:end);

sig_num=size(Y,2);
cnt=1;
while cnt<3
    X=[];
    for i=1:sig_num
        y=Y(:,i);
        x=omp_svd(Dic_est,y);
        X=[X,x];
    end
    
    for i=1:atm_num
        
        x=X(i,:);
        w=find(x);
        omiga=zeros(sig_num,length(w));
        for j=1:length(w)
            
            omiga(w(j),j)=1;
        end
        d=Dic_est(:,i);
        E=Y-Dic_est*X+d*x;
        E_R=E*omiga;
        [U,V,D]=svd(E_R);
        Dic_est(:,i)=U(:,1);
        x_R=V(1,1)*D(:,1);
        
        X(i,w)=x_R;
        
    end
    cnt=cnt+1;
    
end
dic_find=[];



%% 

for i=1:atm_num
d=dic(:,i);
p=abs(d'*Dic_est);
[val,pos]=max(p);
dic_find=[dic_find,Dic_est(:,pos)];
Dic_est(:,pos)=zeros(sig_len,1);
end
rel_vec=diag(dic'*dic_find)

