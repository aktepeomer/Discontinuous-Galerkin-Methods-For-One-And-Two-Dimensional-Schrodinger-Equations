function [V] = b2_Vlinear(example,xpnts,ypnts,nx,ny,basisno,h,plypow)
%genel olarak yanlis duzeltilmesi lzim hem amin de hem burda

ci=complex(0,1);


V=zeros(nx*ny*basisno);
xm=(xpnts(1:nx)+xpnts(2:nx+1))/2;
ym=(ypnts(1:nx)+ypnts(2:nx+1))/2;
for i=1:nx
    for j=1:ny
        c=(i-1)*ny+j; % order of subsquares
        for k=1:basisno %k-th test fnc
            for l=1:basisno %l-th test fnc
                f=@(x,y) (sin(x).^2.*sin(y).^2).*... %potential
                    ((x-xm(i))./(h/2)).^(plypow(k,1)+plypow(l,1)).*...
                    ((y-ym(j))./(h/2)).^(plypow(k,2)+plypow(l,2));
                
                V((c-1)*basisno+k,(c-1)*basisno+l)...
                    =integral2(f,xpnts(i),xpnts(i+1),ypnts(j),ypnts(j+1));
            end
        end
    end
end




end

