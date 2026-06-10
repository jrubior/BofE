function y = LOmeToSigmaQ(x,info)


n=info.nvar;
n2 = n*n;


D=h_tilde(reshape(diag(x(1:n)),n,n),info);

% index_all  = 1:n2;
% index_diag = linspace(1,n2,n);
% index_nondiag = setdiff(index_all,index_diag);

vecL0= x(1+n:n2+n);
L0 = reshape(vecL0,n,n);
Sigma = L0*D*(L0*D)';

Q = inv(chol(Sigma)')*L0*D;





y = [vec(Sigma);vec(Q)];

