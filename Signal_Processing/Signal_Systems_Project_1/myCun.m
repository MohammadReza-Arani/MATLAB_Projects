% A GENERALAZED CONVOLUTION COMPUTING CODE IN MATLAB WITHOUT USING MATLAB BUILTIN FUNCTION conv(x,h)
function Y=myCun(p,d)
%x=input('Enter x:   ')
%h=input('Enter h:   ')
m=length(p);
n=length(d);
X=[p,zeros(1,n)]; 
H=[d,zeros(1,m)]; 
for i=1:n+m-1
    Y(i)=0;
for j=1:m
if(i-j+1>0)
Y(i)=Y(i)+X(j)*H(i-j+1);
else
end
end
end
%stem(Y);
%ylabel('Y[n]');
%xlabel('----->n');
%title('Convolution of Two Signals without conv function');
