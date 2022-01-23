%Çó¸ùx^5-1=0;

function y=rootnewton(z)
if (z==0)
    y=0;
    return;
end
for l=1:1:100000
    y=z-(z^5-1)/(5*z^4);
    if(abs(y-z)<1.0e-3)
    break;
    end
    z=y;
end