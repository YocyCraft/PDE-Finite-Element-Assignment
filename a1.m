function a=a1(i,j,n)%一次有限元方程系数计算
a=0;
ix=mod(i-1,n)+1;
iy=ceil(i/n);
jx=mod(j-1,n)+1;
jy=ceil(j/n);
if(ix==jx && iy==jy)
    a=4;
elseif(abs(jx-ix)==1 && iy==jy)||(abs(jy-iy)==1 && ix==jx)
    a=-1;
end
end