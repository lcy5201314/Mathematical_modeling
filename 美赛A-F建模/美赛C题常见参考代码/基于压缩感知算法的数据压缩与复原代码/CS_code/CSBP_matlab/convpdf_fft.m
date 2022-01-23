function c=convpdf_fft(a,b,epsilon)

tmp=fft((b));
ind=find(abs(tmp)<epsilon);
tmp(ind)=epsilon;
c=a.*tmp;
