function [soln, RHS]=c_InitialConditionConvert(k,nx,ny,xx,yy,h,Alocal,plypow,example)
%in descriptions when we say "x^n" we usually mean {(x-xmid)/(h/2)}^n
%these are test and base functions
%for xmid=(x_{i+1}+x_{i})/2 and h=x_{i+1}-x_{i}
%x_{i}'s are points we take on x-axis. Same for y


basisno=(k+2)*(k+1)/2; % total number of basis polys with degree k or less

rhsx=zeros(nx,k+1); rhsy=zeros(ny,k+1); % \int_{subsquare}(x^n * exact)
xm=(xx(1:nx) +xx(2:nx+1))/2; 	% for each subsquare mid of x 
ym=(yy(1:ny) +yy(2:ny+1))/2; 	% for each subsquare mid of y 
xl=xx(1:nx);     % for each subsquare the left x point
xr=xx(2:nx+1);   % for each subsquare the right x point
yl=yy(1:ny);     % for each subsquare the left y point
yr=yy(2:ny+1);	% for each subsquare the right y point
ci=complex(0,1);% I probably used i before so "i" need this
RHS=zeros(nx*ny*basisno,1);%rhs after integrating by test functions
soln=zeros(nx*ny*basisno,1);%coefficients of galerkin polynomials

if example==0%0000000000000000000000000000000000000000000000000000000000000    
    for n=0:k %n for x^n, k fo P^k poly spc
        rhsx(:,n+1) = ci *(h/2)^(-n) *(ci)^(-n) *exp(-ci*xm) .*...
            ( igamma(n+1,ci*(xr-xm))-igamma(n+1,ci*(xl-xm)) );
        rhsy(:,n+1) = ci *(h/2)^(-n) *(ci)^(-n) *exp(-ci*ym) .*...
            ( igamma(n+1,ci*(yr-ym))-igamma(n+1,ci*(yl-ym)) );
    end
    for i=1:nx
        for j=1:ny
            c=(i-1)*ny+j; %order of subsquares
            for n=1:basisno
                RHS((c-1)*basisno+n)=rhsx(i,plypow(n,1)+1)*rhsy(j,plypow(n,2)+1);
            end
            %compute each region separately they are separate anyway
            %this is much faster that having a big matrix
            soln((c-1)*basisno+1:c*basisno)=...
                linsolve(Alocal,RHS((c-1)*basisno+1:c*basisno));
        end
    end
elseif example==1%111111111111111111111111111111111111111111111111111111111
    for i=1:nx
        for j=1:ny
            c=(i-1)*ny+j; %order of subsquares
            for n=1:basisno
                f=@(x,y) 2^0.5 *exp( ci*(x+y)).*...
                      ((x-xm(i))./(h/2)).^(plypow(n,1)).*...
                      ((y-ym(j))./(h/2)).^(plypow(n,2));
                RHS((c-1)*basisno+n)=integral2(f,xx(i),xx(i+1),yy(j),yy(j+1));
%                 RHS((c-1)*basisno+n)=integral2(f,xx(i),xx(i+1),yy(j),yy(j+1),'AbsTol',1e-16);
            end
            soln((c-1)*basisno+1:c*basisno)=...
                linsolve(Alocal,RHS((c-1)*basisno+1:c*basisno));
        end
    end
elseif example==2%222222222222222222222222222222222222222222222222222222222
    for i=1:nx
        for j=1:ny
            c=(i-1)*ny+j; %order of subsquares
            for n=1:basisno
                f=@(x,y) sin(x).*sin(y).*...
                      ((x-xm(i))./(h/2)).^plypow(n,1).*...
                      ((y-ym(j))./(h/2)).^plypow(n,2);
                RHS((c-1)*basisno+n)=integral2(f,xx(i),xx(i+1),yy(j),yy(j+1));
                %RHS((c-1)*basisno+n)=integral2(f,xx(i),xx(i+1),yy(j),yy(j+1),'AbsTol',1e-12);
            end
            soln((c-1)*basisno+1:c*basisno)=...
                linsolve(Alocal,RHS((c-1)*basisno+1:c*basisno));
        end
    end
end

end

