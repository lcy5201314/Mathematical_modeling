function g=Gold(n,length)
poly_x=zeros(1,26);
poly_x(1)=1;
poly_x(26)=1;
poly_x(23)=1;
poly_y=zeros(1,26);
poly_y(1)=1;
poly_y(23:26)=1;

hpn1=commsrc.pn;
hpn1.GenPoly=poly_x;
hpn1.InitialStates(1:24)=bitget(n,1:24);

hpn2=commsrc.pn;
hpn2.GenPoly=poly_y;
hpn2.InitialStates(1:24)=ones(1,24);

hpn1.NumBitsOut=length;
hpn2.NumBitsOut=length;

pn1=generate(hpn1);
pn2=generate(hpn2);
g1=xor(pn1,pn2);
pn1=gf(pn1,2);
pn2=gf(pn2,2);
poly11=gf([1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0],2);
poly12=gf([1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 ],2);
pn12=conv(pn1,poly11);
pn22=conv(pn2,poly12);
g2=pn12(1:length)+pn22(1:length);
g=(1-2*g1)+j*(1-2*logical(g2.x));