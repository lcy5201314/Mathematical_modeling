function [ mn,mm,mdata ] = myMin( a )
%����������Сֵ�������Ӧ�±�
mdata=inf;
[n,m] = size(a);
for i=1:n
    for j=1:m
        if(a(i,j)<mdata)
            mdata = a(i,j);
            mn = i;
            mm = j;
        end
    end
end
end

