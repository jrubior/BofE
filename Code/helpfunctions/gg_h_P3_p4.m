function y = gg_h_P3_p4(x, info)
%
% Mapping from (B, Sigma, Q) to (B, Sigma, P3, p4)
%

nvar=info.nvar;
npredetermined=info.npredetermined;

B=reshape(x(1:nvar*npredetermined),npredetermined,nvar);
Ssigma=reshape(x(nvar*npredetermined+1:nvar*(npredetermined+nvar)),nvar,nvar);
Qdraw = reshape(x(end-nvar*nvar+1:end),nvar,nvar);

nex = info.nex;
nlag = info.nlag;


tmp=zeros(nvar,nvar);
for l=1:nlag
        tmp = tmp + B(nex+(l-1)*nvar+1:nex+l*nvar,:)';
end

temp_xBS =(eye(nvar)-tmp)\(chol(Ssigma)');

xBSt     = temp_xBS(1,:)./norm(temp_xBS(1,:));
e = eye(nvar);
y = xBSt'-e(:,nvar);
H = eye(nvar)-2*(y*y')/(y'*y);
P = H\Qdraw; 

P3 = P(1:3,1:3);
p4 = P(4,4);

y=[x(1:nvar*(npredetermined+nvar)); vec(P3);p4];