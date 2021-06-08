%一次Lagrange元求解偏微分方程
n=floor(1/input('请输入网格尺寸'));
h=1/n;%网格尺寸
g=h/2;
n=n-1;
A=zeros(n^2,n^2);
B=zeros(n^2,1);
s_mu=h^2/2;%三角形元的面积
for i=1:n^2
    for j=[i,i+1,i+n]
        if(j<=n^2)
            A(i,j)=a1(i,j,n);%计算方程系数
        end
    end
end
for i=2:n^2
    for j=[i-n,i-1]
        if(j>=1)
            A(i,j)=A(j,i);
        end
    end
end
for j=1:n^2
    B(j,1)=b1(j,h,n);%计算右端项系数
end
u=A\B;
uh=zeros(n+2,n+2);%将列向量u重新排布为插值近似uh
for i=1:n+2
    for j=1:n+2
        if(i==1||i==n+2||j==1||j==n+2)
            uh(i,j)=0;
        else
            uh(i,j)=u((j-2)*n+i-1);
        end
    end
end
err0=0;%计算L2(H0)度量下的误差
for i=1:n+1
    for j=1:n+1
        x=i*h;
        y=j*h;
        err0=err0+(u0(x-g,y)-(uh(i+1,j+1)+uh(i,j+1))/2)^2+(u0(x,y-g)-(uh(i+1,j+1)+uh(i+1,j))/2)^2+2*(u0(x-g,y-g)-(uh(i+1,j+1)+uh(i,j))/2)^2+(u0(x-g,y-h)-(uh(i+1,j)+uh(i,j))/2)^2+(u0(x-h,y-g)-(uh(i,j+1)+uh(i,j))/2)^2;
    end
end
err0=sqrt(err0*s_mu/3);
err1=0;%计算H1度量下的误差
for i=1:n+1
    for j=1:n+1
        x=i*h;
        y=j*h;
        uh_grad_left=[uh(i+1,j)-uh(i,j),uh(i,j+1)-uh(i,j)]/h;
        uh_grad_right=[uh(i+1,j+1)-uh(i,j+1),uh(i+1,j+1)-uh(i+1,j)]/h;
        err1=err1+sqnorm(u0_grad(x-h,y-g)-uh_grad_left)+sqnorm(u0_grad(x-g,y-h)-uh_grad_left)+sqnorm(u0_grad(x-g,y-g)-uh_grad_left)+sqnorm(u0_grad(x,y-g)-uh_grad_right)+sqnorm(u0_grad(x-g,y)-uh_grad_right)+sqnorm(u0_grad(x-g,y-g)-uh_grad_right);
    end
end
err1=sqrt(err1*s_mu/3);
disp(err0);
disp(err1);