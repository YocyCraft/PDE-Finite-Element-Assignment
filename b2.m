function b=b2(j,h,n)
jx=mod(j-1,n)+1;
jy=ceil(j/n);
x=jx*h/2;
y=jy*h/2;
s_mu=h^2/2;
%���ֵ���ֵ����
if ~(mod(jx,2)==0 && mod(jy,2)==0)%��i���е�
    b=s_mu*f(x,y)*2/3;
else
    b=0;
end
