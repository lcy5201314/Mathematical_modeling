function p2=ifftshift_shri(p1)

p2=zeros(size(p1));
l=length(p1);
s=(l-1)/2;
p2(s+2:l)=p1(1:s);
p2(1:s+1)=p1(s+1:l);
