function spreadingBySir()
    A=load('test.in');
    % node number
    N=size(A,1);   
    %��Ⱦ����
    irate=0.3;
    %�ָ�����
    rrate=1;
    %��ʼʱ�ڵ��״̬��,��ʼʱֻ�нڵ�1Ϊ��Ⱦ״̬�������Ķ�Ϊ�׸�Ⱦ״̬  
    start_node=1;
    %����ͼ�Ĺ�����ȵ�ԭ�������в�������
    BFSspreading(A,N,start_node,irate,rrate);
end

function BFSspreading(A,N,start_node,irate,rrate)
%����ͷ
head=1;            
%����β����ʼ����Ϊ�գ�tail==head
tail=1;            
%��ͷ�м����ȾԴ�ڵ�
queue(head)=start_node;      
%������չ
head=head+1;  

%��Ⱦ�ڵ��б� 
infection=start_node;  
%�ָ��ڵ��б�  
recover=[];
%�׸�Ⱦ�ڵ��б�
for i=1:N
    %��ʼʱ��start_nodeΪ��Ⱦ״̬
    if i==start_node
        %-1��ʾ�ýڵ��Ѿ����б���ɾ��
        susceptible(i)=-1;
    end
    %��ʼʱ������start_nodeΪ��Ⱦ״̬�⣬�����ڵ㶼�����׸�Ⱦ״̬
     susceptible(i)=i;
end

%��ʼ���չ����������˳�����ھӽڵ㴫��
%�ж϶����Ƿ�Ϊ��
while tail~=head   
    %ȡ��β�ڵ� 
    i=queue(tail);  
    %����ýڵ㲻���Ƴ��б�֮��
    if isempty(find(recover==i,1))
            for j=1:N
             %����ڵ�j�뵱ǰ�ڵ�i�������ҽڵ�j���ڸ�Ⱦ�б���
            if A(i,j)==1 && isempty(find(infection==j,1))   
                 infection_random=rand(1);
                 if infection_random < irate
                    %�½ڵ�����
                    queue(head)=j;  
                    %��չ����
                    head=head+1;   
                    %���½ڵ�j�����Ⱦ�б�
                    infection=[infection j]; 
                    
                    %���׸�Ⱦ�ڵ��б���ɾ���ýڵ�,����Ϊ-1
                    [row,col,v] = find(susceptible==j) ;
                    susceptible(col)=-1;
                    susceptible(find(susceptible==-1))=[];                    
                 end
            end
        end
        %����Ⱦ�Ľڵ㰴���ʼ���ָ��ڵ��б�  
        recover_random=rand(1);
        if infection_random < rrate
            %�ָ�
            recover=[recover i];  
            %�Ӹ�Ⱦ�б���ɾ��
            [row,col,v] = find(infection==i) ;
            infection(col)=-1;
            infection(find(infection==-1))=[];
        end
        tail=tail+1; 
        
    end %end if  isempty(find(recover==i,1)
end %end while

%�ֱ���ʵ���ڵ��״̬
infection
susceptible
recover
end

