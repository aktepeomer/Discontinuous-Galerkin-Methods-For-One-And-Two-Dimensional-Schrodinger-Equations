function [WW,QQ,A,plypow] = b1_IntegralMatrices(k,nx,ny,beta,h,dt)
%This part returns first 3 integrals matrices
%in descriptions when we say x we actually mean (x-xmid)/(h/2) 
%for xmid=(x_{i+1}+x_{i})/2 and h=x_{i+1}-x_{i}
%x_{i}'s are grid points we take on x-axis. Same for y 

basisno=(k+2)*(k+1)/2; % total number of basis polys

plypow=zeros(basisno,2);% powers of x, powers of y
ply_xpow=zeros(basisno,3);  ply_ypow=zeros(basisno,3);
% powers of x, powers of y, coeff coming from derivation
count=1;
for p=0:k %total power of the polinom
    for pp=0:p %power of y
        plypow  (count,:)=[p-pp  ,pp  ];
        ply_xpow(count,:)=[p-pp-1,pp  ,(p-pp)];%some values are -1. it is useful ...
        ply_ypow(count,:)=[p-pp  ,pp-1,(  pp)];%when p_{x}(x,y)=0 and p_{y}(x,y)=0
        count=count+1;
    end
end

%Main matrices (integral) in the paper with this order
A =zeros(basisno);B1=zeros(basisno);B2=zeros(basisno);
C =cell(2,4,3); %1-2(in-ot); 1:4(sides); 1:3(3 integrals)
% basis poly is indepentent from the region so compute just for one region
for i=1:basisno %i is for test function, gives rows
    for j=1:basisno %j is for soln. u, gives columns
        %for dt term integrals
        kk=plypow(i,1)+plypow(j,1);%total power of x
        ll=plypow(i,2)+plypow(j,2);%total power of y
        A(i,j)= (h/2)^2 *(1-(-1)^(kk+1)) *(1-(-1)^(ll+1)) /((kk+1)*(ll+1));
        
        % for the integral with nabla psi \cdot nabla basis
        % d_x terms
        if ply_xpow(i,1)==-1 || ply_xpow(j,1)==-1 % derivative is 0-> integral is 0
            B1(i,j)=0;
        else
            k_x=ply_xpow(i,1)+ply_xpow(j,1);%total power of x for dx poly
            l_x=ply_xpow(i,2)+ply_xpow(j,2);%total power of y for dx poly
            coef=ply_xpow(i,3)*ply_xpow(j,3);
            B1(i,j)= coef *(1-(-1)^(k_x+1)) *(1-(-1)^(l_x+1)) /((k_x+1)*(l_x+1));
        end
        % d_y terms
        if ply_ypow(i,2)==-1 || ply_ypow(j,2)==-1 % derivative is 0-> integral is 0
            B2(i,j)=0;
        else
            k_y=ply_ypow(i,1)+ply_ypow(j,1);%total power of x for dy poly
            l_y=ply_ypow(i,2)+ply_ypow(j,2);%total power of y for dy poly
            coef=ply_ypow(i,3)*ply_ypow(j,3);
            B2(i,j)= coef *(1-(-1)^(k_y+1)) *(1-(-1)^(l_y+1)) /((k_y+1)*(l_y+1));
        end
        
        
        %line integrals:
        %first inside integrals: right,top,left,bottom 1,2,3,4
        kk=plypow(i,1)+plypow(j,1);%total power of x
        ll=plypow(i,2)+plypow(j,2);%total power of y
        C{1,1,1}(i,j)= (beta/h) *( 1)^kk *(h/2) *(1-(-1)^(ll+1)) / (ll+1);
        C{1,3,1}(i,j)= (beta/h) *(-1)^kk *(h/2) *(1-(-1)^(ll+1)) / (ll+1);
        C{1,2,1}(i,j)= (beta/h) *( 1)^ll *(h/2) *(1-(-1)^(kk+1)) / (kk+1);
        C{1,4,1}(i,j)= (beta/h) *(-1)^ll *(h/2) *(1-(-1)^(kk+1)) / (kk+1);
        %first outside integrals: right,top,left,bottom 1,2,3,4
        C{2,1,1}(i,j)=-(beta/h) *(-1)^plypow(j,1) *(h/2) *(1-(-1)^(ll+1)) / (ll+1);
        C{2,3,1}(i,j)=-(beta/h) *(-1)^plypow(i,1) *(h/2) *(1-(-1)^(ll+1)) / (ll+1);
        C{2,2,1}(i,j)=-(beta/h) *(-1)^plypow(j,2) *(h/2) *(1-(-1)^(kk+1)) / (kk+1);
        C{2,4,1}(i,j)=-(beta/h) *(-1)^plypow(i,2) *(h/2) *(1-(-1)^(kk+1)) / (kk+1);
        
        
        %second inside/outside integrals: right,top,left,bottom 1,2,3,4
        %right and left
        if ply_xpow(j,1)==-1 % derivative is 0-> integral is 0
            C{1,1,2}(i,j)=0;    C{1,3,2}(i,j)=0;
            C{2,1,2}(i,j)=0;    C{2,3,2}(i,j)=0;
        else
            k_x=plypow(i,1)+ply_xpow(j,1);%total power of x
            l_x=plypow(i,2)+ply_xpow(j,2);%total power of y
            coef_x=ply_xpow(j,3);
            %inside
            C{1,1,2}(i,j)=-0.5 *coef_x *( 1)^k_x *(1-(-1)^(l_x+1)) /(l_x+1);
            C{1,3,2}(i,j)= 0.5 *coef_x *(-1)^k_x *(1-(-1)^(l_x+1)) /(l_x+1);
            %outside
            C{2,1,2}(i,j)=-0.5 *coef_x *(-1)^ply_xpow(j,1) *(1-(-1)^(l_x+1)) /(l_x+1);
            C{2,3,2}(i,j)= 0.5 *coef_x *(-1)^plypow(i,1)   *(1-(-1)^(l_x+1)) /(l_x+1);
        end
        %top and bottom
        if ply_ypow(j,2)==-1 % derivative is 0-> integral is 0
            C{1,2,2}(i,j)=0;    C{1,4,2}(i,j)=0;
            C{2,2,2}(i,j)=0;    C{2,4,2}(i,j)=0;
        else
            k_y=plypow(i,1)+ply_ypow(j,1);%total power of x
            l_y=plypow(i,2)+ply_ypow(j,2);%total power of y
            coef_y=ply_ypow(j,3);
            %inside
            C{1,2,2}(i,j)=-0.5 *coef_y *( 1)^l_y *(1-(-1)^(k_y+1)) /(k_y+1);
            C{1,4,2}(i,j)= 0.5 *coef_y* (-1)^l_y *(1-(-1)^(k_y+1)) /(k_y+1);
            %outside
            C{2,2,2}(i,j)=-0.5 *coef_y *(-1)^ply_ypow(j,2) *(1-(-1)^(k_y+1)) /(k_y+1);
            C{2,4,2}(i,j)= 0.5 *coef_y* (-1)^plypow(i,2)   *(1-(-1)^(k_y+1)) /(k_y+1);
        end
        
        %third inside/outside integrals: right,top,left,bottom 1,2,3,4
        %right and left
        if ply_xpow(i,1)==-1 % derivative is 0-> integral is 0
            C{1,1,3}(i,j)=0;    C{1,3,3}(i,j)=0;
            C{2,1,3}(i,j)=0;    C{2,3,3}(i,j)=0;
        else
            k_x=ply_xpow(i,1)+plypow(j,1);%total power of x
            l_x=ply_xpow(i,2)+plypow(j,2);%total power of y
            coef_x=ply_xpow(i,3);
            %inside
            C{1,1,3}(i,j)=-0.5 *coef_x *( 1)^k_x *(1-(-1)^(l_x+1)) /(l_x+1);
            C{1,3,3}(i,j)= 0.5 *coef_x *(-1)^k_x *(1-(-1)^(l_x+1)) /(l_x+1);
            %outside
            C{2,1,3}(i,j)= 0.5 *coef_x *(-1)^plypow(j,1)   *(1-(-1)^(l_x+1)) /(l_x+1);
            C{2,3,3}(i,j)=-0.5 *coef_x *(-1)^ply_xpow(i,1) *(1-(-1)^(l_x+1)) /(l_x+1);
        end
        %top and bottom
        if ply_ypow(i,2)==-1 % derivative is 0-> integral is 0
            C{1,2,3}(i,j)=0;    C{1,4,3}(i,j)=0;
            C{2,2,3}(i,j)=0;    C{2,4,3}(i,j)=0;
        else
            k_y=ply_ypow(i,1)+plypow(j,1);%total power of x
            l_y=ply_ypow(i,2)+plypow(j,2);%total power of y
            coef_y=ply_ypow(i,3);
            %inside
            C{1,2,3}(i,j)=-0.5 *coef_y *( 1)^l_y *(1-(-1)^(k_y+1)) /(k_y+1);
            C{1,4,3}(i,j)= 0.5 *coef_y* (-1)^l_y *(1-(-1)^(k_y+1)) /(k_y+1);
            %outside
            C{2,2,3}(i,j)= 0.5 *coef_y *(-1)^plypow(j,2)   *(1-(-1)^(k_y+1)) /(k_y+1);
            C{2,4,3}(i,j)=-0.5 *coef_y* (-1)^ply_ypow(i,2) *(1-(-1)^(k_y+1)) /(k_y+1);
        end
        
    end
end

Amat=zeros(nx*ny*basisno);Bmat=zeros(nx*ny*basisno);Cmat=zeros(nx*ny*basisno);
for i=1:nx
    for j=1:ny
        c=(i-1)*ny+j; % order of subsquares 
        ct=mod(c+ny-1,nx*ny)+1;% top-bottom-left-right c neighbours
        cb=mod(c-ny-1,nx*ny)+1;
        cl=(i-1)*ny+mod(j-1-1,ny)+1;
        cr=(i-1)*ny+mod(j+1-1,ny)+1;
        
        Amat((c-1)*basisno+1:(c)*basisno, (c-1)*basisno+1:(c)*basisno)...
            = A;
        Bmat((c-1)*basisno+1:(c)*basisno, (c-1)*basisno+1:(c)*basisno)...
            = B1+B2;
        Cmat((c-1)*basisno+1:(c)*basisno, (c-1)*basisno+1:(c)*basisno)...
            = Cmat((c-1)*basisno+1:(c)*basisno, (c-1)*basisno+1:(c)*basisno)...
             +C{1,1,1}+C{1,2,1}+C{1,3,1}+C{1,4,1}...
             +C{1,1,2}+C{1,2,2}+C{1,3,2}+C{1,4,2}...
             +C{1,1,3}+C{1,2,3}+C{1,3,3}+C{1,4,3};
        %in the order right top left bottom 
        Cmat((c-1)*basisno+1:(c)*basisno, (cr-1)*basisno+1:(cr)*basisno)...
      = Cmat((c-1)*basisno+1:(c)*basisno, (cr-1)*basisno+1:(cr)*basisno)...
            +C{2,1,1}+C{2,1,2}+C{2,1,3};
        Cmat((c-1)*basisno+1:(c)*basisno, (ct-1)*basisno+1:(ct)*basisno)...
      = Cmat((c-1)*basisno+1:(c)*basisno, (ct-1)*basisno+1:(ct)*basisno)...
            +C{2,2,1}+C{2,2,2}+C{2,2,3};
        Cmat((c-1)*basisno+1:(c)*basisno, (cl-1)*basisno+1:(cl)*basisno)...
      = Cmat((c-1)*basisno+1:(c)*basisno, (cl-1)*basisno+1:(cl)*basisno)...
            +C{2,3,1}+C{2,3,2}+C{2,3,3};
        Cmat((c-1)*basisno+1:(c)*basisno, (cb-1)*basisno+1:(cb)*basisno)...
      = Cmat((c-1)*basisno+1:(c)*basisno, (cb-1)*basisno+1:(cb)*basisno)...
            +C{2,4,1}+C{2,4,2}+C{2,4,3};
    end
end


ci=complex(0,1);
WW=ci*Amat*(2/dt) -Bmat-Cmat;
QQ=ci*Amat*(2/dt);


end

