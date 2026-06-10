function y = gg_h_P3_p4_inv(x, info)
%
% Mapping from (B, Sigma, P3,p4) to (B, Sigma, Q)
%

nvar=info.nvar;
npredetermined=info.npredetermined;

break0=nvar*npredetermined;
B=reshape(x(1:break0),npredetermined,nvar);
break1=break0+nvar*nvar;
Ssigma=reshape(x(break0+1:break1),nvar,nvar);
break2=break1+(nvar-1)*(nvar-1);
P3=reshape(x(break1+1:break2),nvar-1,nvar-1);
p4=reshape(x(end),1,1);

nex = info.nex;
nlag = info.nlag;

P = [P3 zeros(3,1)
     zeros(1,3)  p4];


tmp=zeros(nvar,nvar);
for l=1:nlag
        tmp = tmp + B(nex+(l-1)*nvar+1:nex+l*nvar,:)';
end

temp_xBS =(eye(nvar)-tmp)\(chol(Ssigma)');
xBSt     = temp_xBS(1,:)./norm(temp_xBS(1,:));

e = eye(nvar);
y = xBSt'-e(:,nvar);
H = eye(nvar)-2*(y*y')/(y'*y);
Qdraw = H*P; 

q=reshape(Qdraw,nvar*nvar,1);

y=[x(1:nvar*(npredetermined+nvar)); q];