function y = LOmeToBSigmaQ(x,info)


n=info.nvar;
n2 = n*n;


%D=reshape(diag(x(1:n)),n,n);
D=chol(reshape(diag(x(1:n)),n,n));

%index_all  = 1:n2;
% index_diag = linspace(1,n2,n);
% index_nondiag = setdiff(index_all,index_diag);

%L0 = eye(n);
%vecL0 = L0(:);
%vecL0(index_nondiag) = x(1+n:n2);
vecL0= x(1+n:n2+n);
L0 = reshape(vecL0,n,n);


Sigma = L0*D*(L0*D)';

Q = inv(chol(Sigma)')*L0*D;

L1p = nan(n*info.nlag,n);

if info.nex==1
for ell=1:info.nlag
    %L1p(1+(ell-1)*info.nvar:ell*info.nvar,:)= reshape(x(ell*n2+1:(ell+1)*n2),n,n);
    L1p(1+(ell-1)*info.nvar:ell*info.nvar,:)= reshape(x(n+ell*n2+1:n+(ell+1)*n2),n,n);
end

c   = x(end-n+1:end,1)';
else
    
for ell=1:info.nlag
    %L1p(1+(ell-1)*info.nvar:ell*info.nvar,:)= reshape(x(ell*n2+1:(ell+1)*n2),n,n);
    L1p(1+(ell-1)*info.nvar:ell*info.nvar,:)= reshape(x(n+ell*n2+1:n+(ell+1)*n2),n,n);
end
end

B = nan(info.m,n);




for ell=1:info.nlag


    tmp = zeros(n);

    for k=1:(ell-1)

        tmp = tmp + (L1p(1+((ell-k)-1)*info.nvar:(ell-k)*info.nvar,:)/L0)'*B(info.nex+1+(k-1)*info.nvar:info.nex+k*info.nvar,:);
    
    end
    B(info.nex+1+(ell-1)*info.nvar:info.nex+ell*info.nvar,:) = (L1p(1+(ell-1)*info.nvar:ell*info.nvar,:)/L0)'-tmp;

end

if info.nex==1

B(1,:) = c*D*Q'*chol(Sigma);

end


y = [vec(B);vec(Sigma);vec(Q)];