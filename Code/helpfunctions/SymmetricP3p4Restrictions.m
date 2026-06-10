function y = SymmetricP3p4Restrictions(x, info)

nvar=info.nvar;
npredetermined=info.npredetermined;

y=zeros(nvar*(nvar-1)/2 + (nvar-1)*((nvar-1)+1)/2 + 1,1);

% Sigma must be symmetric
Sigma = reshape(x(nvar*npredetermined+1:nvar*npredetermined+nvar*nvar),nvar,nvar);
k=1;
for i=1:nvar
    for j=i+1:nvar
        y(k)=Sigma(i,j) - Sigma(j,i);
        k=k+1;
    end
end

% P3 must be orthogonal
Q = reshape(x(nvar*(npredetermined+nvar)+1:end-1),nvar-1,nvar-1);
for i=1:(nvar-1)
    y(k)=Q(:,i)'*Q(:,i) - 1.0;
    k=k+1;
    for j=i+1:(nvar-1)
        y(k)=Q(:,i)'*Q(:,j);
        k=k+1;
    end
end

% % p4 must be of norm one
p4 = x(end);
y(k)= norm(p4)^2 - 1.0;


