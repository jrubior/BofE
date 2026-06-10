function y = BSigmaQToOBJ(x,info)


nvar=info.nvar;
m=info.m;
nvar2=info.nvar^2;

B=reshape(x(1:m*nvar),m,nvar);
Sigma=reshape(x(m*nvar+1:(m+nvar)*nvar),nvar,nvar);
Q=reshape(x((m+nvar)*nvar+1:end),nvar,nvar);

D = diag(diag(chol(Sigma)'*Q));

L0= chol(Sigma)'*(Q); %/D

p = (m-info.nex)/nvar;

A     = cell(p,1);
A0    = inv(L0)';
Aplus = B*A0;

L=cell(p,1);

for i=1:p
    A{i}=Aplus(info.nex+(i-1)*nvar+1:info.nex+i*nvar,:);
    X = A0\A{i};
    for j=1:i-1
        X = X + L{j}*A{i-j};
    end
    L{i}=X/A0;
end

if info.nex==1
y = nan(nvar2+nvar2*info.nlag+nvar,1);
else
y = nan(nvar2+nvar2*info.nlag,1);
end

break0=nvar;
y(1:break0)=reshape(diag(D),nvar,1);

L0norm = (inv(A0)')/D;
vecL0norm  = L0norm(:);

index_diag = linspace(1,nvar2,nvar);

vecL0norm(index_diag) = [];

break1             = break0 + length(vecL0norm);
y(break0+1:break1) = vecL0norm;


for i=1:p
    y(break1+1+(i-1)*nvar2:break1+i*nvar2)=reshape((L{i}')/D,nvar2,1);
end

if info.nex==1
break2=break1+p*nvar2;
break3=break2 + nvar;
y(break2+1:break3)=reshape(Aplus(1,:),nvar,1);

end





vecL0 = L0(:);

index_diag = linspace(1,nvar2,nvar);

vecL0(index_diag) = [];

y(nvar+1:nvar2)=vecL0;





if info.nex==1
c = inv(D)'*Q'*inv(chol(Sigma))'*(B(1,:)');
y((info.nlag+1)*nvar2+1:end)=c;

end

