function a=a2(i,j,n)%��������Ԫ����ϵ������
a=0;
ix=mod(i-1,n)+1;
iy=ceil(i/n);
jx=mod(j-1,n)+1;
jy=ceil(j/n);
if(mod(ix,2)==0 && mod(iy,2)==0)%��i�Ƕ���
    if (ix==jx && iy==jy)%i��j��ͬһ��
        a=4;
    elseif((ix==jx && abs(iy-jy)==1)||(iy==jy && abs(ix-jx)==1))%i��j�����ڵ�
        a=-4/3;
    elseif((ix==jx && abs(iy-jy)==2)||(iy==jy && abs(ix-jx)==2))%i��j���һ��
        a=1/3;
    end
else%��i���е�
    if (ix==jx && iy==jy)%i��j��ͬһ��
        a=16/3;
    elseif(ix==jx && abs(iy-jy)==1)||(iy==jy && abs(ix-jx)==1)%i��j�����ڵ�
        a=-4/3;
    end
end