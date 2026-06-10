function y = pphi_P3_p4(x,info)

% Mapping from (L0,L1,L2,L3,L_inf,c) to (B,Sigma,P3,p4).  If (L0,L1,L2,L3,L_inf,c) satisfies the zero
% restrictions, then pphi_P3_p4_inv(pphi_P3_p4(L0,L1,L2,L3,L_inf,c)) = (L0,L1,L2,L3,L_inf,c)


nex = info.nex;
nlag = info.nlag;
nvar=info.nvar;
%npredetermined=info.npredetermined;

break0 = nvar*nvar;
L0    = reshape(x(1:break0,1),nvar,nvar);
break1 = break0+nvar*nvar;
L1    = reshape(x(1+break0:break1,1),nvar,nvar);
break2 = break1+nvar*nvar;
L2    = reshape(x(1+break1:break2,1),nvar,nvar);
break3 = break2+nvar*nvar;
L3    = reshape(x(1+break2:break3,1),nvar,nvar);
break4 = break3+nvar*nvar;
L_inf    = reshape(x(1+break3:break4,1),nvar,nvar);
break5 = break4+nvar;
c      = reshape(x(1+break4:break5,1),1,nvar);

A0         = inv(L0)';
Ssigmadraw = (A0*A0')\eye(nvar);

A1 = (L1/L0)'*A0;
A2 = (L2/L0)'*A0 - (L1/L0)'*A1;
A3 = (L3/L0)'*A0 - (L2/L0)'*A1 - (L1/L0)'*A2;

%Linf = (A0'-sumA1A4t)\eye(nvar);
% inv(Linf) = (A0'-sumA1A4t);
% sumA1A4t = A0'-inv(Linf);
% A4t = A0'-inv(Linf)-sumA1A3';
% A4  = (A0'-inv(Linf)-sumA1A3')';
A4 = (A0'-inv(L_inf)-(A1'+A2'+A3'))';

Aplus = [c;A1;A2;A3;A4];

Bdraw      = Aplus/A0;

tmp=zeros(nvar,nvar);
for l=1:nlag
        tmp = tmp + Bdraw(nex+(l-1)*nvar+1:nex+l*nvar,:)';
end

temp_xBS =(eye(nvar)-tmp)\(chol(Ssigmadraw)');
Qdraw = temp_xBS\L_inf;

xBSt     = temp_xBS(1,:)./norm(temp_xBS(1,:));
e = eye(nvar);
y = xBSt'-e(:,nvar);
H = eye(nvar)-2*(y*y')/(y'*y);
P = H\Qdraw; 

P3 = P(1:3,1:3);
p4 = P(4,4);




 y=[vec(Bdraw);vec(Ssigmadraw); vec(P3);p4];

