%%�������˿ڵĻ����Ⱦ�������˿���������ɵĴ�Ⱦ
%a=input('ÿ������ÿ����Ч�Ӵ�����ƽ��ֵ�� ');
%b=input('ÿ�챻�������˱����� ');
%c=input('ÿ��ÿ�������˿�ת�Ʊ���: ');
%y0=input('��ʼ��Ⱦ�߱����� ');
%w=input('��������кţ� ');
disp('��������֮������ڹ�ϵ�� ');
A=ans
N=size(A,2);
a=1;
b=0.3;
c=0.004;
y0=0.02;
w=N;
i=zeros(1,N);
i(w)=y0;
for t=1:16    %����
    for j=1:N  %9���������˿ڵĻ����Ⱦ
        i(j)=((tanh((a/2 - b/2)*(1 + (2*atanh((2*a*i(j))/(a - b) - 1))/(a - b))) + 1)*(a - b))/(2*a);
    end
    for m=1:N-1  %���ڳ��м��˿�����
        for n=(m+1):N
            if 1==A(m,n)
                i(m)=i(m)-c*i(m)+c*i(n);
                i(n)=i(n)-c*i(n)+c*i(m);
            end
        end
    end
end
disp('16���������еĸ�Ⱦ��Ϊ��');
i

%%�����������֮���·������
D=A;
D(find(D==0))=inf;
for j=1:N
    D(j,j)=0;
end
for k=1:N
    for m=1:N
        for j=1:N
            if(D(m,j)>D(m,k)+D(k,j))
                D(m,j)=D(m,k)+D(k,j);
            end
        end
    end
end
a=find(D==inf);
D(a)=0;
mm=max(D(w,:));  %mmΪ�������������w���������е���Զ·��ֵ
b=find(D==0);
D(b)=mm;
for j=1:N
    D(j,j)=0;
end
disp('��������֮���·������Ϊ: ');
D

%%����Ⱦ�����������ϴ�Ⱦ��Χ�ĺ���ͼ
ii=zeros(1,mm); %%ii: ���и�Ⱦ��������뱬����·�����ȵĹ�ϵ
for j=1:mm
    flag=0; 
    for m=1:N-1
        if j==D(N,m)
            ii(j)=ii(j)+i(m);
            flag=flag+1;
        end
    end
    ii(j)=ii(j)/flag;
end
ii
bar(ii);
title('����SIS�Ĵ�Ⱦ�����������ϴ�Ⱦ��Χ������ͼ');
xlabel('��Ⱦ·������');
ylabel('��Ⱦ����');


