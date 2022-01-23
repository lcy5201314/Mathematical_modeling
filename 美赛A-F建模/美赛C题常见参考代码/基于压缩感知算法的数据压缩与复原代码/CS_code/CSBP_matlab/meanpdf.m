function [m]=meanpdf(pdf, xx);

pdf=pdf/sum(pdf);
m=sum (xx.*pdf);
