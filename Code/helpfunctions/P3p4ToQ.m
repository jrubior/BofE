function q = P3p4ToQ(x, info, Bdraw, Ssigmadraw)
%Z_IRF, nvar, W)
%
%  w=[w(1); ... ; w(n)] where the dimension of w(j) is n-(j-1+z(j)) > 0. 
%
%  Z_IRF{j} - z(j) x n matrix of full row rank. 
%      Z_IRF{j} = Z{j} * F(f_h^{-1}(B,Sigma,I_n)).
%
%  W{j} - (n-(j-1+z(j))) x n matrix.
%
%  q=[q(1); ... ; q(n)] where the dimension of q(j) is n.
%
%  norm(q(j)) == 1 if and only if norm(x(j)) == 1.
%

nex = info.nex;
nlag = info.nlag;
nvar=info.nvar;


P3 = reshape(x(1:(nvar-1)*(nvar-1),1),nvar-1,nvar-1);
p4 = x(end,1);

P = [P3 zeros(3,1)
     zeros(1,3)  p4];


tmp=zeros(nvar,nvar);
for l=1:nlag
        tmp = tmp + Bdraw(nex+(l-1)*nvar+1:nex+l*nvar,:)';
end

temp_xBS =(eye(nvar)-tmp)\(chol(Ssigmadraw)');
xBSt     = temp_xBS(1,:)./norm(temp_xBS(1,:));

e = eye(nvar);
y = xBSt'-e(:,nvar);
H = eye(nvar)-2*(y*y')/(y'*y);
Qdraw = H*P; 


q=reshape(Qdraw,nvar*nvar,1);
