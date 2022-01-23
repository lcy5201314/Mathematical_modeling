%function [p,A]=randomgraph()
function suijitu()
disp('���ͼ���ɲ���1,2,3��4')
disp('1��ʾ �����������ȽϷ����ܹ����ɵı���ΪN*(N-1)/2*alph');
disp('2��ʾ �������򷨣��ܹ����ɵı���ΪN*(N-1)/2*alph,����һ���Ľ�С�ĸ��ʶԱ����������');
disp('3��ʾ �����������Ƚϣ�����Ҫ���ܹ��ı���ΪN*(N-1)/2*alph');
disp('4��ʾ ���ַ����ܹ����ɵı���ΪN*(N-1)/2*alph');
pp=input('���������ͼ���ɲ���1,2,3��4:');
% N=input('����ͼ�нڵ������ĿN��');
% alph=input('����ͼ�бߵ�ƽ�����Ӷ�alph:  ');
% beta=input('�����ߵ�ƽ�����ȵĲ���beta:  ');
N=100;
alph=0.25;
beta=0.3;
randData=rand(2,N)*1000;
x=randData(1,:);
y=randData(2,:);
p=lianjiegailv(x,y,alph,beta,N);
switch  pp
    case 1
        A=bian_lianjie1(p,N,alph);
    case 2
         relink=input('��������������ӵĸ���:');
         A=bian_lianjie2(p,N,alph,relink);
    case 3
         A=bian_lianjie3(p,N,alph);
    case 4
         A=bian_lianjie4(p,N,alph);
    otherwise
         disp('The number pp you input is wrong');
         return;
end

plot(x,y,'r.','Markersize',18);
hold on;
for i=1:N 
    for j=i+1:N
        if A(i,j)~=0
            plot([x(i),x(j)],[y(i),y(j)],'linewidth',1);
            hold on;
        end
    end
end
hold off
 
[C,aver_C]=Clustering_Coefficient(A);
[DeD,aver_DeD]=Degree_Distribution(A);
[D,aver_D]=Aver_Path_Length(A);   
 disp(['�����ͼ�ı���Ϊ��',int2str(sum(sum(A))/2)]); 
 disp(['�����ͼ��ƽ��·������Ϊ��',num2str(aver_D)]);  %%������������������
 disp(['�����ͼ�ľ���ϵ��Ϊ��',num2str(aver_C)]);
 disp(['�����ͼ��ƽ����Ϊ��',num2str(aver_DeD)]);   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ú���������ڵ����ӱߵĸ���
function p=lianjiegailv(x,y,alph,beta,N)
d=zeros(N);
for i=1:N
    for j=1:N
        d(i,j)=sqrt((x(i)-x(j))^2+((y(i)-y(j)))^2);
    end
end
L=max(max(d));
for i=1:N
    for j=1:N
        p(i,j)=alph*exp(-d(i,j)/beta/L);
    end
        p(i,i)=0;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���ɻ���1����[0,1]�ھ���������Ƚϣ���p(i,j)>random_data,�����ӽڵ�i,j.
% ֱ���ܹ����ɵı���ΪN*(N-1)/2*alph
function A=bian_lianjie1(p,N,alph)   %  ����ֵDΪ�ڽӾ���
A=zeros(N);num=0;
for k=1:inf
    for i=1:N
        for j=1:N
            random_data=rand(1,1);
            if p(i,j)>=random_data&A(i,j)==0
               A(i,j)=1;A(j,i)=1;
                num=num+1;
                if num>=N*(N-1)/2*alph
                   return ;
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ɻ���2�������ʴӴ�С�������Ӹ�������ǰ��Ľڵ�ԣ�ֱ���ܹ����ɵı���ΪN*(N-1)/2*alph
%% ��һ���Ľ�С����������������ʵ��һ���̶ȵ��������   �����⣡��������������
function A=bian_lianjie2(p,N,alph,relink)
A=zeros(N);
p1=reshape(tril(p),[1,N*N]);
[p2,px]=sort(p1,'descend');
M=ceil(N*(N-1)/2*alph)
for k=1:M
    [m,n]=ind2sub(size(p),px(k));  %���±�������Ϊ˫�±�����
    A(m,n)=1;A(n,m)=1;
end
[m,n]=find(tril(A));               %��һ���ĸ������������
for i=1:length(m)
    p1=rand(1,1);
    if relink>p1                          %����������������ʴ������ɵ�������������������
         A(m(i),n(i))=0;A(n(i),m(i))=0;   %�ȶϿ�ԭ���ıߣ������ѡ��һ������֮����  
         A(m(i),m(i))=inf;
         n1=find(A(m(i),:)==0);      
         random_data=randint(1,1,[1,length(n1)]);
         nn=n1(random_data);
         A(m(i),nn)=1;A(nn,m(i))=1;
         A(m(i),m(i))=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ɻ���3����[0,1]�ھ���������Ƚϣ���p(i,j)>random_data,�����ӽڵ�i,j��
%% ����Ҫ���ܹ��ı���ΪN*(N-1)/2*alph
function A=bian_lianjie3(p,N,alph);     
A=zeros(N);
for i=1:N
    for j=1:N
         random_data=rand(1,1);
         if p(i,j)>=random_data&A(i,j)==0
            A(i,j)=1; A(j,i)=1;
         end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ɻ���4�������ʹ�һ�������ö��ַ�ѡ�����ӵıߣ�ֱ�����ɱ���ΪN*(N-1)/2*alph��
function A=bian_lianjie4(p,N,alph)
A=zeros(N);
p1=reshape(p,1,N*N)./sum(sum(p));   
pp=cumsum(p1);%���ۼƸ���
k=0;
while  k<N*(N-1)/2*alph           %���ö��ַ�ѡ��һ��������
     random_data=rand(1,1);
     aa=find(pp>=random_data);jj=aa(1); % �ڵ�jj��Ϊ�ö��ַ�ѡ��Ľڵ�
     j=ceil(jj/N);i=jj-(j-1)*N;             %�ѵ��±�������Ϊ˫�±������������ú���ind2sub(siz,IND)
    % [i,j=ind2sub(size(p),jj);
    if A(i,j)==0
        A(i,j)=1;A(j,i)=1;
        k=k+1;
    end   
end

% function A=bian_lianjie4(p,N,alph) 
% %%���ɻ���4������˫�ζ��ַ���ÿ��ѡ��һ�������ӵĽڵ㣬�ٽ������ڵ�������ֱ�����ɱ���N*(N-1)/2*alph��
% %  �������������ѡ��һ�У����ö��ַ�ѡ������ӱߵ�һ���ڵ�i�����ڵ�i������ö��ַ�ѡ����һ���ڵ�j�������������ڵ��γ�һ����
% A=zeros(N);
% for i=1:N
%     p(i,:)=p(i,:)./sum(p(i,:));          %��ÿһ�еĸ��ʹ�һ��
% end
% k=0;
% pp=cumsum(p,2);                          %�Թ�һ�����ÿһ�еĸ������ۼӺ�
% while k<N*(N-1)/2*alph
%     kk=randint(1,1,[1,N]);              %���ѡ��һ��
%     random_data1=rand(1,1);
%     ii=find(pp(kk,:)>=random_data1);
%     i=ii(1);                            %���ö��ַ�ѡ������ӱߵ�һ���ڵ�i
%     random_data2=rand(1,1);
%     jj=find(pp(i,:)>=random_data2);
%     j=jj(1);                            %���ڵ�i������ö��ַ�ѡ����һ���ڵ�j
%     if A(i,j)==0
%       A(i,j)=1;A(j,i)=1;
%       k=k+1;
%     end
% end
    






