u=@(t) double(t >= 0);
n = -5:1:5;
h= sinc(2*(pi)*n).*(u(n+4)-u(n-5));
x=u(n)-u(n-2);

tic
y1=conv(h,x);
t1=toc
figure
stem(convindicesme(n,n),y1,'linewidth',2)
tic
y2=myCun(h,x);
t2=toc
figure
stem(convindicesme(n,n),y2,'linewidth',2)



