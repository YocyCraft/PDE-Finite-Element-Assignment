function b=b1(j,h,n)
jx=mod(j-1,n)+1;
jy=ceil(j/n);
x=jx*h;
y=jy*h;
g=h/2;
s_mu=h^2/2;
b=s_mu*(f(x+g,y)+f(x,y+g)+f(x-g,y+g)+f(x-g,y)+f(x,y-g)+f(x+g,y-g))/3;%积分的数值计算

%b=h^2*f(x,y);