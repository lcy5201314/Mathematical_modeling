function p2=shiftpdf(p1, s, delta, mo)

p1=p1(length(p1):-1:1);
p1=[zeros(1, mo), p1, zeros(1,mo)];
sind=round(-s/delta);
sind=sind+(length(p1)+1)/2;
lind=sind-(mo-1)/2;
rind=sind+(mo-1)/2;
p2=p1(lind:rind);
p2=p2/sum(p2);
