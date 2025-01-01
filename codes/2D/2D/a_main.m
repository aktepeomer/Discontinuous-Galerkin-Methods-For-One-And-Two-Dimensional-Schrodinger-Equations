function [uffff] = a_main(nx,k)
tic
example=1;
%0:linear,1:potential(u)=2|u|^2, 2:potential(u)=2|u|^2-sin^2(x)*sin^2(y)u
%domains are different be careful
beta=0;
tfinal=1;   dt=0.002; %final time, dt

k=2; %k is P^k(x,y) poly spc
nx=10*3;


ny=nx;%number of points at x-axis is nx+1

%boundary points of x and y axis
if example==0 ||example==1
    xa=0;   xb=2*pi;
    ya=0;   yb=2*pi;
elseif example==2
    xa=0;   xb=pi;
    ya=0;   yb=pi;
end
ci=complex(0,1);
dx=(xb-xa)/nx;  dy=(yb-ya)/ny;
x=xa:dx:xb;     y=ya:dy:yb;     %points taken for [x_0,x_1]...[x_n-1,x_n]
t=0:dt:tfinal;  nt=round(tfinal/dt)+1; %for time; nt=final time index
h=dx;

basisno=(k+2)*(k+1)/2;

soln=zeros(nx*ny*basisno,nt);
% % soln=zeros(nx*ny*basisno,1);

[WW,QQ,Alocal,plypow]=b1_IntegralMatrices(k,nx,ny,beta,h,dt);
[soln(:,1),RHS]=c_InitialConditionConvert(k,nx,ny,x,y,h,Alocal,plypow,example);

if example==2
    EE=b2_Vlinear(example,x,y,nx,ny,basisno,h,plypow);
    QQ=QQ+EE;
    clear EE
end

M=linsolve( WW, QQ );
%size and condition numbers
% disp('size and condition number')
% size(WW)
% cond(WW)
% condest(WW)
% 1/(rcond(WW))
clear WW QQ
M=2*M-eye(nx*ny*basisno);

if example==0%0000000000000000000000000000000000000000000000000000000000000
    for n=1:nt-1
        soln(:,n+1)=M*soln(:,n);
% %         soln(:)=M*soln(:);
    end
elseif example==1 || example==2%121212121212121212121212121212121212121212
    for n=1:nt-1
        ustar=b3_VnonlinearStrang(Alocal,nx,ny,soln(:,n),plypow,basisno,dt,x,y);
% %         ustar=b3_VnonlinearStrang(Alocal,nx,ny,soln(:),plypow,basisno,dt,x,y);
        u2str=M*ustar;
        soln(:,n+1)=b3_VnonlinearStrang(Alocal,nx,ny,u2str,plypow,basisno,dt,x,y);
% %         soln(:,1)=b3_VnonlinearStrang(Alocal,nx,ny,u2str,plypow,basisno,dt,x,y);
    end
end
clear M

totalerr=zeros(nt,1);

% error=zeros(nx*ny,nt);
% error=zeros(nt,1);

for time=nt:nt
    err=zeros(nx,ny);
    for i=1:nx
        for j=1:ny
            c=(i-1)*ny+j;% order of subsquares
            aa=repmat( soln((c-1)*basisno+1:c*basisno,time), 1, basisno);
% %             aa=repmat( soln((c-1)*basisno+1:c*basisno), 1, basisno);
            apap=sum( (aa.*aa').*Alocal,'all');
            
            if example==0%0000000000000000000000000000000000000000000000000
                exex=h*h;
                exap=exp(-ci*2*t(time))*...
                    soln((c-1)*basisno+1:c*basisno,time)' * RHS((c-1)*basisno+1:c*basisno);
% %                 exap=exp(-ci*2*t(time))*...
% %                     soln((c-1)*basisno+1:c*basisno)' * RHS((c-1)*basisno+1:c*basisno);
                
            elseif example==1%111111111111111111111111111111111111111111111
                exex=2*h*h;
                exap=exp( ci*2*t(time))*...
                    soln((c-1)*basisno+1:c*basisno,time)' * RHS((c-1)*basisno+1:c*basisno);
% %                 exap=exp( ci*2*t(time))*...
% %                     soln((c-1)*basisno+1:c*basisno)' * RHS((c-1)*basisno+1:c*basisno);
                
            elseif example==2%222222222222222222222222222222222222222222222
                exex=0.5*(x(i+1)-sin(x(i+1))*cos(x(i+1))-x(i)+sin(x(i))*cos(x(i)))*...
                     0.5*(y(j+1)-sin(y(j+1))*cos(y(j+1))-y(j)+sin(y(j))*cos(y(j)));
                exap=exp(-ci*2*t(time))*...
                    soln((c-1)*basisno+1:c*basisno,time)' * RHS((c-1)*basisno+1:c*basisno);
            end
            err(i,j)=exex-2*real(exap)+apap;
%             [exex,real(exap),apap]
        end
    end
    totalerr(time)=(abs(sum(err,'all')))^0.5;
%     totalerr(time)
%     error(:,time)=err(:);
%     error(time)=err(1,1);
end
% totalerr
% totalerr(1:nt).'
% error.'
totalerr(nt)
toc

% save('k4icin.mat','soln')

uffff=totalerr(end);
end

