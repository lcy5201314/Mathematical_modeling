function A=qujian(B) %%%%%��ÿһ��ҩƷ������������

[m,n]=size(B);
for i=1:m
    k=B(i,4); %%%��
    c=B(i,2); %%%��
    g=B(i,3); %%%��
    dkg=sqrt(k^2+g^2); %%%��߶Խ���
    dck=sqrt(c^2+k^2); %%%���Խ���
    A(i,1)=B(i,1); %%���
    A(i,2)=k+2;%%����
    A(i,3)=min([2*k,dkg,dck]);%%����
end
    