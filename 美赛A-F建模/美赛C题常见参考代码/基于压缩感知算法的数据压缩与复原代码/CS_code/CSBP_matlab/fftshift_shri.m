function p2=fftshift_shri(p1)

p2=zeros(size(p1));
l=length(p1);
s=(l-1)/2;
p2(s+1:l)=p1(1:s+1);
p2(1:s)=p1(s+2:l);
