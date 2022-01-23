function []=Animation1()
n = 20;
s = .01;
x = rand(n,1)-0.5;
y = rand(n,1)-0.5;
h = plot(x,y,'.');
axis([-1 1 -1 1]);
set(h,'EraseMode','xor','MarkerSize',18)
while 1
   drawnow;
   x = x + s*randn(n,1);
   y = y + s*randn(n,1);

   set(h,'XData',x,'YData',y)
end