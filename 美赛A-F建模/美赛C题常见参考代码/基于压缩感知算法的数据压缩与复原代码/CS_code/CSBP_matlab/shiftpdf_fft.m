function p4=shiftpdf_fft(p1, shiftamt, delta, mo)

l=length(p1);
p3=zeros(1,l);
p4=p3;
s=(l-1)/2;
p3(s+1:l) = p1(1:s+1);
p3(1:s) = p1(s+2:l);
p3=p3(l:-1:1);
p3b=zeros(1,mo*2+l);
p3b(mo+1:mo+l)=p3;
sind=round(-shiftamt/delta);
sind=max(min(sind,mo),-mo); % solves out of bounds crash, dB, 11.22.2008
sind=sind+s+1+mo;
p2=p3b(sind-s:sind+s);
p2=p2/sum(p2);
p4(s+2:l)=p2(1:s);
p4(1:s+1)=p2(s+1:l);
