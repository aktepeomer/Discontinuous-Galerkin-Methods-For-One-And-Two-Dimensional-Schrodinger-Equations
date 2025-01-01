function [unew] = explicit_part(FD2,FD4,uold,k,N)
% returns RK methods explicit part
% FD2,FD4 are cells with dim=k+1,1
% each cell is a matrix
% each matrix obtained by (i-1)-th poly. with coeff |a_i|^2(or 4)
% i-th cell=> l,m entry is "int of p_j p_m (p_i)^2(or 4)"

cnj_uold=conj(uold);
FD=zeros(N*(k+1));


for i=0:k
    for j=0:k
        dumvec1=zeros(N,1);
        for subints=1:N
            dumvec1(subints)= uold( (subints-1)*(k+1)+i+1 )...
                *cnj_uold( (subints-1)*(k+1)+j+1 );
        end
        coefof2=repelem(dumvec1,(k+1),N*(k+1));
        FD=FD+FD2{i+1,j+1}.*coefof2;
        
        for ii=0:k
            for jj=0:k
                dumvec2=zeros(N,1);
                for subints=1:N
                    dumvec2(subints)= uold( (subints-1)*(k+1)+ii+1 )...
                        *cnj_uold( (subints-1)*(k+1)+jj+1 );
                end
                coefof4=coefof2.*repelem(dumvec2,(k+1),N*(k+1));
                FD=FD+FD4{i+1,j+1,ii+1,jj+1}.*coefof4;
            end
        end
    end
end

unew=FD*uold;

end