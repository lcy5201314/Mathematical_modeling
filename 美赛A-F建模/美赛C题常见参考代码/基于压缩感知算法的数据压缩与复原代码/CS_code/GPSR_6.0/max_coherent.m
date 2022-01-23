function [u,max] = max_coherent( A )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[m,n]=size(A);
max=0;
for i=1:n
    for j=1:n
       v= A(:,i).*A(:,j);
       if v>max
           max=v;
           vi=norm(A(:,i),2);
              vj=norm(A(:,i),2);
       end
    end

end
u=max/vi*vj;

