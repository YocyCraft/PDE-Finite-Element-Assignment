function a=a2(i,j,n)%二次有限元方程系数计算
a=0;
ix=mod(i-1,n)+1;
iy=ceil(i/n);
jx=mod(j-1,n)+1;
jy=ceil(j/n);
if(mod(ix,2)==0 && mod(iy,2)==0)%点i是顶点
    if (ix==jx && iy==jy)%i与j是同一点
        a=4;
    elseif((ix==jx && abs(iy-jy)==1)||(iy==jy && abs(ix-jx)==1))%i与j是相邻点
        a=-4/3;
    elseif((ix==jx && abs(iy-jy)==2)||(iy==jy && abs(ix-jx)==2))%i与j相隔一点
        a=1/3;
    end
else%点i是中点
    if (ix==jx && iy==jy)%i与j是同一点
        a=16/3;
    elseif(ix==jx && abs(iy-jy)==1)||(iy==jy && abs(ix-jx)==1)%i与j是相邻点
        a=-4/3;
    end
end