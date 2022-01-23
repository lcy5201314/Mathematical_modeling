function k = nz( x ,delta)
%This function is find the sparisty of signal x
%   Detailed explanation goes here

n=size(x);
for i=1:n
    if abs(x(i))<delta
        x(i)=0;
    else 
       x(i)=1;
       
    end
end
k=sum(x);
fprintf(1,'the length x is %10.3d\n',n);
fprintf(1,'the sparisty of x is %10.3d\n',k);

