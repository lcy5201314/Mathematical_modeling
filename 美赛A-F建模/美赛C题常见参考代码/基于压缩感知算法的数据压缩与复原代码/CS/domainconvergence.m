clc;
clear;
a=[1,0,0,0,0,-1];x=roots(a)


 A=x(1);


 B=x(2);

 C=x(3);

 D=x(4);

 E=x(5);

 h=0.01;

for a=-2:h:2
       for b=-2:h:2
         z=a+b*i;   
         y=rootnewton(z);
         if (abs(y-A)<1.0e-3)
             plot(a,b,'r');
%              hold on
%          elseif (abs(y-B)<1.0e-3)
%              plot(a,b,'y');
%              hold on
%          elseif (abs(y-C)<1.0e-3)
%              plot(a,b,'g');
%              hold on    
%           elseif (abs(y-D)<1.0e-3)
%              plot(a,b,'b');
%              hold on
%           elseif (abs(y-E)<1.0e-3)
%              plot(a,b,'p');
%              hold on     
         end
       end
end      
