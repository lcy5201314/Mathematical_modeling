function [m,sig,mp]=meanvarmaxpdf(pdf, xx);

pdf=pdf/sum(pdf);
[mm,ind]=max(pdf);
mp=xx(ind);
m=sum (xx.*pdf);
e2=sum(xx.*xx.*pdf);
v2=e2-m*m;
sig=sqrt(v2);
