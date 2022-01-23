function dispvec_anderrors (v, phi, phisign, measvec, dispind, x);

l=length(dispind);
s='[';
for i=1:l
  s=sprintf('%s %7.2f',s,v(dispind(i)));
end
s=sprintf('%s]',s);
mv=encoder(phi, phisign, v);
er=norm(measvec-mv);
s=sprintf('%s (%7.2f) (%7.2f)',s, er, norm(x-v));
disp(s);
