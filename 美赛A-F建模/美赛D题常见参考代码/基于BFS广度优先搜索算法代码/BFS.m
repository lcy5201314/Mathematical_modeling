%%��������
% zhilu=[
%     1 2   ;
%     1 6   ;
%     1 7   ;
%     2 3   ;
%     2 8   ;
%     3 4   ;
%     3 9   ;
%     4 5   ;
%     4 10  ;
%     5 6   ;
%     5 11  ;
%     6 12  ;
%     7 8   ;
%     7 12  ;
%     8 9   ;
%     9 10  ;
%     10 11 ;
%     11 12 ; 
%    ];
zhilu1=[ 
    0 3 ;
    1 2 ;
    2 3 ;
    2 4 ;
    3 5 ;
    5 7 ;
    5 9 ;
    7 6 ;
    9 8 
    ];

%%������֧·����ת��Ϊ�ڽӾ���
[m1,n1]=size(zhilu1);
zhilu=zhilu1+ones(m1,n1);
n=max(max(zhilu(:,1:2)));                 %��ȡ֧·�ڵ���
G=zeros(n);       
for i=1:m1
  m2=zhilu(i,1);
  n2=zhilu(i,2);
  G(m2,n2)=1;
  G(n2,m2)=1;
end
%%Ѱ�����һ������������Ķ���
W=zeros(1,n);                            %�����ź�Ľڵ㣬�ڵ�˳���С��������
l=0;
v=1;
a1=find(G(v,:)==1);                      %Ѱ�����һ������������ڵ㲢���
G(v,a1)=2;                               
G(a1,v)=2;
W(a1)=l+1;
S1=union(v,a1);
l=l+1;
%%Ѱ������Ϊl�Ķ����������δ����ŵĶ��㼯��
while ~isempty(G==1)

    a1=find(G(S1,:)==1);
    t=length(S1);
    d=[];
    for i=1:length(a1)
        if a1(i)/t>floor(a1(i)/t)
            t2=floor(a1(i)/t)+1;
        else
            t2=floor(a1(i)/t);
        end                              %col
        if isempty(intersect(d,t2))
            d=union(d,t2);
        end
    end
    d1= setdiff(d,S1);
    %���ҵ��Ķ��㼯�Ͻ��б��
    if isempty(d1)
        break;
    else 
        W(d1)=l+1;
        G1=G(S1,:);
        G1(a1)=2;
        G(S1,:)=G1;
        G(:,S1)=G1';
        S1=union(S1,d1);
        l=l+1;
    end
    
end




















