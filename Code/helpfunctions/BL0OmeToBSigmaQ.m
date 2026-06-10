function y = BL0OmeToBSigmaQ(x,info)


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

vecB= x(n2+n+1:end,1);

y = [vecB;vec(Sigma);vec(Q)];