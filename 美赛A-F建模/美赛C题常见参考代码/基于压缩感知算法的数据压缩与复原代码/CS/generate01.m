function y=genarate01(x);
x=randn(2,4);
if length(size(x))<2 
    error('the size of x must less than 3')
end

y=zeros(size(2,4));
for i=1:2
    for j=1:4
        if abs(x(i,j))<1
            y(i,j)=0;
        else
            y(i,j)=1;
        end
    end
end

            
        

