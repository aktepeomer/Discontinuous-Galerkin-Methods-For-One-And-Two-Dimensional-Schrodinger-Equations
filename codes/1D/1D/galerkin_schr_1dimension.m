function [yenierrl2] = galerkin_schr_1dimension(N,k)
%we consider equidistant intervals for now
N=10*1;     % number of intervals
k=2;        % max poly degree
dt=1/100;  % time steps

a=-pi;
b= pi;


alpha1=0;   %flux paremeters
alpha2=0;
beta1=0;
beta2=0;


x_bndry=(a:(b-a)/N:b);                        %N+1 many points
% x_center=(x_bndry(1:N)+x_bndry(2:N+1))/2;   %N many points
h=(b-a)/N;


dummymat= repmat(0:k,k+1,1) + repmat( (0:k).',1,k+1);
dummymat2=( ones(k+1) + repmat(-1,k+1).^dummymat )/2;
dummymat1=dummymat+ones(k+1);
dummymat4=dummymat-ones(k+1);
if k>0
    dummymat4(1:2,:)=1;
elseif k==0
    dummymat4(1:1,:)=1;
end
dummymat3=repmat( ((0:k).*(-1:k-1)).' ,1,k+1);

% explaniation is on the paper
A = h     ./ dummymat1 .* dummymat2;
A = A*complex(0,1);
B = (4/h) ./ dummymat4 .* dummymat3 .* dummymat2;
FA=kron( eye(N), A );
Ainv=A^(-1);
% FB=FA^(-1)*kron( eye(N), B );
FB=kron( eye(N), Ainv*B );



% nu parts of C part:
% at x_j+1/2 and x_j-1/2
% rows are same
% poly value at x_j-1/2. rows are same:
C_nu_plus= repmat(-1,k+1).^repmat( (0:k).',1,k+1);
% poly value at x_j+1/2. rows are same. all positive:
C_nu_minus= ones(k+1);
% derivative of poly at x_j+1/2. rows are same. all positive:
C_nu_x_minus= (2/h)*repmat((0:k).',1,k+1);
% derivative of poly at x_j-1/2. rows are same:
C_nu_x_plus= C_nu_x_minus.*C_nu_plus*-1;


% for each "constant" term there are 4 matrices:
% p_x^- *uhat at x_j+1/2 coeff.
AAA=(0.5+alpha2)*C_nu_x_minus.*C_nu_plus.';
AAB=(0.5-alpha2)*C_nu_x_minus.*C_nu_minus.';
AAC=beta2*       C_nu_x_minus.*C_nu_x_plus.';
AAD=beta2*       C_nu_x_minus.*C_nu_x_minus.';
% p_x^+ *uhat at x_j-1/2 coeff.
ABA=(0.5+alpha2)*C_nu_x_plus.*C_nu_plus.';
ABB=(0.5-alpha2)*C_nu_x_plus.*C_nu_minus.';
ABC=beta2*       C_nu_x_plus.*C_nu_x_plus.';
ABD=beta2*       C_nu_x_plus.*C_nu_x_minus.';
% p_x^- *utilt at x_j+1/2 coeff.
BAA=(0.5+alpha1)*C_nu_minus.*C_nu_x_plus.';
BAB=(0.5-alpha1)*C_nu_minus.*C_nu_x_minus.';
BAC=beta1*       C_nu_minus.*C_nu_plus.';
BAD=beta1*       C_nu_minus.*C_nu_minus.';
% p_x^+ *utilt at x_j-1/2 coeff.
BBA=(0.5+alpha1)*C_nu_plus.*C_nu_x_plus.';
BBB=(0.5-alpha1)*C_nu_plus.*C_nu_x_minus.';
BBC=beta1*       C_nu_plus.*C_nu_plus.';
BBD=beta1*       C_nu_plus.*C_nu_minus.';

FC_diog=kron( eye(N),-(AAB-AAD)+(ABA+ABC)+(BAB-BAD)-(BBA+BBC));
upp_diog=full(gallery('tridiag',N,0,0,1));
dow_diog=full(gallery('tridiag',N,1,0,0));
FC_upp=kron(upp_diog,-(AAA+AAC)+(BAA+BAC));
FC_dow=kron(dow_diog,(ABB-ABD)-(BBB+BBD));
FC=FC_diog+FC_upp+FC_dow;
FC(1:k+1, (N-1)*(k+1)+1:N*(k+1))=(ABB-ABD)-(BBB+BBD);
FC((N-1)*(k+1)+1:N*(k+1), 1:k+1)=-(AAA+AAC)+(BAA+BAC);
FC=FA^(-1)*FC;
% for diffusion term. the whole interval
% we include the constant terms too.
% ??? is this correct?
% because of the form on the RK papere
FBC=-(FB+FC);



% for |u|^2 + |u|^4 = f(|u|^2)
% these have to be multiplied by |a_i||a_j|(and 4)
D2=cell(k+1,k+1);
D4=cell(k+1,k+1,k+1,k+1);
% for whole intvl.s.must be multiplied by different a_i'son each subint
FD2=cell(k+1,k+1);
FD4=cell(k+1,k+1,k+1,k+1);
for i=0:k
    for j=0:k
        if mod(i+j,2)==0
            somemat=dummymat2;
        else
            somemat=ones(k+1)-dummymat2;
        end
        D2{i+1,j+1}=h./(dummymat1+ i*ones(k+1)+j*ones(k+1)).*somemat;
        D2{i+1,j+1}=-1*Ainv*D2{i+1,j+1};
        FD2{i+1,j+1}=kron( eye(N), D2{i+1,j+1} );
        for ii=0:k+1
            for jj=0:k+1
                if mod(i+ii+j+jj,2)==0
                    somemat=dummymat2;
                else
                    somemat=ones(k+1)-dummymat2;
                end
                D4{i+1,j+1,ii+1,jj+1}=h./(dummymat1+ ...
                    i*ones(k+1)+j*ones(k+1)+...
                    ii*ones(k+1)+jj*ones(k+1)).*somemat;
                D4{i+1,j+1,ii+1,jj+1}=-1*Ainv*D4{i+1,j+1,ii+1,jj+1};
                FD4{i+1,j+1,ii+1,jj+1}=...
                    kron( eye(N), D4{i+1,j+1,ii+1,jj+1} );
            end
        end
    end
end
% for sub intervals. need to be multiplied by a_i's
% inverse of A is already applied and send to otherside






% % initial condition. may need to be changed:
% 
% % we add k-1 many points in each interval
% % is this correct???
% x_subs=( a:(b-a)/(N*k):b ).';
% % original initial condition
% IC=exp(complex(0,1)*x_subs);
% % new initial condition
% new_IC=zeros(N*(k+1),1);
% % dd is just y=x (poly.) will define rest of the poly.s by this
% % actually (2x-a-b)/(b-a) but we write here in a simple way
% dd=(-1:2/k:1).';
% % coulumn i is x^i (PS:(2x-a-b)/(b-a) for each interval)
% Mat=zeros(k+1);
% Mat(:,1)=1;
% for m=2:k+1
%     Mat(:,m)=dd.^(m-1);
% end
% invMat=Mat^(-1);
% % IC= Mat * new_IC
% for l=1:N
%     new_IC( (l-1)*(k+1)+1 : l*(k+1) )=...
%         invMat*IC( (l-1)*k+1 : l*k+1 );
% end



% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bu kisim deneme
rhs=zeros(N*(k+1),1);
new_IC=zeros(N*(k+1),1);
ci=complex(0,1);
invA=(-ci*A)^(-1);
for i=1:N
    xm=(x_bndry(i+1)+x_bndry(i))/2;
    
    for n=0:k
        rhs( ((i-1)*(k+1)+1+n) )=exp(ci*xm) *(2/h)^(n)*(-ci)^(-n+1)...
            *( igamma(n+1,ci*(-h)/2)-igamma(n+1,ci*(h)/2) );
    end

    new_IC( (i-1)*(k+1)+1 : i*(k+1) )=invA*rhs( (i-1)*(k+1)+1 : i*(k+1) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bu kisim deneme


% Runge-kutte:

% I use "third order" one in the paper
% is this correct???

gamma=(3+3^0.5)/6;  % given in "third order"
s=2;                % given in "third order"
sigma=s+1;          % true in general

% parameters are given in "third order"
% c=[gamma,1-gamma];
b=[0.5,0.5];       % implicit
bh=[0,0.5,0.5];    % explicit
a=[gamma,0;1-2*gamma,gamma];                       % implicit
ah=[0,0,0 ; gamma,0,0 ; gamma-1,2*(1-gamma),0];    % explicit




phi=zeros((k+1)*N, 1/dt +1 );
phi(:,1)=new_IC;


for tno=2:1/dt +1
    K =zeros((k+1)*N,s);
    Kh=zeros((k+1)*N,sigma);
    phisub=zeros((k+1)*N,s);
    
    Kh(:,1)=explicit_part(FD2,FD4,phi(:,tno-1),k,N);
    for i=1:s
        
        h_sum=zeros((k+1)*N,1);
        n_sum=zeros((k+1)*N,1);
        for j=1:i-1
            h_sum =dt*ah(i+1,j)*Kh(:,j)+h_sum;
            n_sum =dt*a(i,j)   *K(:,j) +n_sum;
        end
            h_sum =dt*ah(i+1,i)*Kh(:,i)+h_sum;
        
        phisub(:,i)=(eye((k+1)*N) - dt*a(i,i)*FBC)^(-1)*...
            (phi(:,tno-1)+h_sum+n_sum);
        
        K(:,i)=FBC*phisub(:,i);
        Kh(:,i+1)=explicit_part(FD2,FD4,phisub(:,i),k,N);
    end
    
    for j=1:s
        phi(:,tno)=phi(:,tno)+dt*b(j)*K(:,j) + dt*bh(j)*Kh(:,j);
    end
    phi(:,tno)=phi(:,tno-1)+phi(:,tno)+dt*bh(sigma)*Kh(:,sigma);
end


t_steps=0:dt:1;
[T,X]=meshgrid(t_steps,x_bndry(1:N));
exact_soln=exp( complex(0,1)* (X+T));


ww=[-1;1];
poly_end=zeros(2,k+1);
for l=0:k
    poly_end(:,l+1)=ww.^l;
end
poly_end_mat=kron(eye(N),poly_end);
int_values=poly_end_mat*phi;

approx=int_values(1:2:2*N,:)+...
    [int_values(2*N,:);int_values(2:2:2*N-2,:)];
approx=approx/2;




% % figure(23)
% % subplot(2,2,1)
% % plot(x_bndry(1:N),real(exact_soln(:,1)),x_bndry(1:N),real(approx(:,1)))
% % subplot(2,2,2)
% % plot(x_bndry(1:N),imag(exact_soln(:,1)),x_bndry(1:N),imag(approx(:,1)))
% % subplot(2,2,3)
% % plot(x_bndry(1:N),real(approx(:,1))-real(exact_soln(:,1)))
% % subplot(2,2,4)
% % plot(x_bndry(1:N),imag(approx(:,1))-imag(exact_soln(:,1)))
% 
% % fff=exact_soln(:,1)-approx(:,1);
% % [h*norm(fff,1),h*norm(fff,2),norm(fff,Inf)]



% % figure(1)
% % surf(T,X,real(approx))
% % figure(2)
% % surf(T,X,imag(approx))
% % 
% % figure(3)
% % surf(T,X,real(exact_soln))
% % figure(4)
% % surf(T,X,imag(exact_soln))
% % 
% % figure(5)
% % surf(T,X,real(approx)-real(exact_soln))
% % figure(6)
% % surf(T,X,imag(approx)-imag(exact_soln))

figure(7)
subplot(2,3,1)
surf(T,X,real(approx))
subplot(2,3,2)
surf(T,X,real(exact_soln))
subplot(2,3,3)
surf(T,X,real(approx)-real(exact_soln))
subplot(2,3,4)
surf(T,X,imag(approx))
subplot(2,3,5)
surf(T,X,imag(exact_soln))
subplot(2,3,6)
surf(T,X,imag(approx)-imag(exact_soln))


fark=approx-exact_soln;
% yenifark=approx(:,end)-exact_soln(:,end)
farkvec=fark(:);

% % figure(12)
% % subplot(2,3,1)
% % subplot(x_bndry(1:N),real(fark(:,1)))
% % subplot(2,3,2)
% % subplot(x_bndry(1:N),real(fark(:,N/2)))
% % subplot(2,3,3)
% % subplot(x_bndry(1:N),real(fark(:,N)))
% % subplot(2,3,4)
% % subplot(x_bndry(1:N),real(fark(:,1)))
% % subplot(2,3,5)
% % subplot(x_bndry(1:N),real(fark(:,N/2)))
% % subplot(2,3,6)
% % subplot(x_bndry(1:N),real(fark(:,N)))

% % 
% % l1norm=h*dt*norm(farkvec,1)
% % l2norm=(h*dt)^(0.5)*norm(farkvec,2)
% % l8norm=max( abs(fark) ,[],'all')


%%%%%yeniiiiiiiiiiii
yenierrl2liste=zeros(N);
yenierrl1liste=zeros(N);
yenierrl8liste=zeros(N);
ci=complex(0,1);
for i=1:N
    exex=h;
    aa=repmat( phi((i-1)*(k+1)+1:i*(k+1),end), 1, (k+1));
    apap=sum( (aa.*aa').*(-ci*A),'all' );
    
    xm=(x_bndry(i+1)+x_bndry(i))/2;
    dumdum=-ci*phi((i-1)*(k+1)+1,end)'*( exp(ci*x_bndry(i+1))-exp(ci*x_bndry(i)) );
    for n=1:k
        dumdum=dumdum+ phi((i-1)*(k+1)+1+n,end)'...
            *exp(ci*xm) *(2/h)^(n)*(-ci)^(-n+1)...
            *( igamma(n+1,ci*(-h)/2)-igamma(n+1,ci*( h)/2) );
    end
    exap=dumdum*exp(ci*t_steps(end));
    yenierrl2liste(i)=abs(exex-2*real(exap)+apap);
    
end

yenierrl2=sum(yenierrl2liste,'all')^0.5
yenierrl1=sum(yenierrl1liste,'all');
yenierrl8=sum(yenierrl8liste,'all');
[yenierrl1,yenierrl2,yenierrl8]

end