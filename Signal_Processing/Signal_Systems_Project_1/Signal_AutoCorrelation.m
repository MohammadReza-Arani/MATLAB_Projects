 clear ; clc
 np=0:5;
 k=ones(size(np));
 ni=convindicesme(np,np);
 
 y=conv(k,k);
 figure
 stem(ni,y,'linewidth',2)
