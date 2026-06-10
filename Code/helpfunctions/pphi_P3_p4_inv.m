function y = pphi_P3_p4_inv(x,info)

% Mapping from (B, Sigma, P3,p4) to (L0,L1,L2,L3,L_inf,c).  Is Sigma is symmetric and
% positive definite, then (L0,L1,L2,L3,L_inf,c) will satisfy the zero restrictions and
% pphi_P3_p4(pphi_P3_p4_inv(B,Sigma,P3,p4)) = (B,Sigma,P3,p4)

nex = info.nex;
nlag = info.nlag;
nvar=info.nvar;
npredetermined=info.npredetermined;

P3 = reshape(x(1+(npredetermined+nvar)*nvar:end-1,1),nvar-1,nvar-1);
p4 = x(end,1);

P = [P3 zeros(3,1)
     zeros(1,3)  p4];

break0     = npredetermined*nvar;
Bdraw      = reshape(x(1:break0,1),npredetermined,nvar);
Ssigmadraw = reshape(x(1+break0:break0+nvar*nvar,1),nvar,nvar);
     
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

L_inf = temp_xBS*Qdraw;

z=f_h_inv([x(1:(npredetermined+nvar)*nvar); vec(Qdraw)],info);

A0=reshape(z(1:nvar*nvar),nvar,nvar);

L0 = inv(A0)';

Aplus=reshape(z(nvar*nvar+1:end),npredetermined,nvar);

c = Aplus(1,:);

A1 = Aplus(2:nvar+1,:);
A2 = Aplus(nvar+2:1+nvar*2,:);
A3 = Aplus(2+nvar*2:1+nvar*3,:);

L1 = (A1/A0)'*L0;
L2 = (A1/A0)'*L1 + (A2/A0)'*L0;
L3 = (A1/A0)'*L2 + (A2/A0)'*L1 + (A3/A0)'*L0;

y=[vec(L0); vec(L1);vec(L2);vec(L3);vec(L_inf);c'];


