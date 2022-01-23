function c=gmean(a,b,gamma, epsilon)

c= (a*(1-gamma)) + (b*(gamma));
c=c/sum(c);
