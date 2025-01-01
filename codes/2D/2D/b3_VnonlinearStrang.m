function [solnnew] = b3_VnonlinearStrang(Alocal,nx,ny,solnold,plypow,basisno,dt,xx,yy)
% plypow=[0,0; 1,0; 0,1];% 2,0; 1,1; 0,2]; 
% basisno=max(size(plypow));
% dt=0.01;
% nx=3;ny=3;  
% xx=0:1:3; yy=xx; 
% solnold=ones(nx*ny*basisno,1)*1;

ci=complex(0,1);
xm=(xx(2:nx+1)+xx(1:nx))./2; 
ym=(yy(2:ny+1)+yy(1:ny))./2;
h=xx(2)-xx(1);


RHS=zeros(nx*ny*basisno,1);%rhs after integrating by test functions
solnnew=zeros(nx*ny*basisno,1);

for i=1:nx
    for j=1:ny
        c=(i-1)*ny+j;% order of subsquares
        
        f=@(x,y) solnold((c-1)*basisno+1).*...
            ((x-xm(i))./(h/2)).^(plypow(1,1)).*...
            ((y-ym(j))./(h/2)).^(plypow(1,2));
        for k=2:basisno
            f=@(x,y) f(x,y) +solnold((c-1)*basisno+k).*...
            ((x-xm(i))./(h/2)).^(plypow(k,1)).*...
            ((y-ym(j))./(h/2)).^(plypow(k,2));
        end
        
        f=@(x,y) f(x,y).*exp( -ci* (-2)*abs((f(x,y))).^2 *(dt/2) );

        for k=1:basisno
            g=@(x,y) f(x,y).*...
                ((x-xm(i))./(h/2)).^(plypow(k,1)).*...
                ((y-ym(j))./(h/2)).^(plypow(k,2));
            RHS((c-1)*basisno+k)=integral2(g,xx(i),xx(i+1),yy(j),yy(j+1));
%             RHS((c-1)*basisno+k)=integral2(g,xx(i),xx(i+1),yy(j),yy(j+1));
        end
        
        solnnew((c-1)*basisno+1:c*basisno)=linsolve(Alocal,RHS((c-1)*basisno+1:c*basisno));
    end
end





end

