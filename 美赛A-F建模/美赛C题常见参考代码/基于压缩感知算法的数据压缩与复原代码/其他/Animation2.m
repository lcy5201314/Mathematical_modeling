function [  ] = Animation2(  )
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here
nframes = 100;
n=20;
s=.02
x = rand(n,1)-0.5;
y = rand(n,1)-0.5;
h = plot(x,y,'yo');
set(h,'MarkerSize',18);
axis([-1 1 -1 1])
for k = 1:nframes
   x = x + s*randn(n,1);
   y = y + s*randn(n,1);
   set(h,'XData',x,'YData',y)
   M(k) = getframe;
end
movie(M,30);