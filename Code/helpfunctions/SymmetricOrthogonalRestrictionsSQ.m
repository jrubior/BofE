function y = SymmetricOrthogonalRestrictionsSQ(x, info)

nvar=info.nvar;
nvar2=nvar*nvar;

y=zeros(nvar*nvar,1);

% Sigma must be symmetric

Sigma=reshape(x(1:nvar2),nvar,nvar);

k=1;
for i=1:nvar
    for j=i+1:nvar
        y(k)=Sigma(i,j) - Sigma(j,i);
        k=k+1;
    end
end

% Q must be orthogonal
Q=reshape(x(nvar2+1:end),nvar,nvar);

for i=1:nvar
    y(k)=Q(:,i)'*Q(:,i) - 1.0;
    k=k+1;
    for j=i+1:nvar
        y(k)=Q(:,i)'*Q(:,j);
        k=k+1;
    end
end