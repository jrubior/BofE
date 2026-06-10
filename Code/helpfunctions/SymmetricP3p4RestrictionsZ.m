function y = SymmetricP3p4RestrictionsZ(x, info)

nvar=info.nvar;
npredetermined=info.npredetermined;

Z=info.Z;
hirfs=info.horizons;

n=size(Z,1);
total_zeros=0;
for j=1:n
    total_zeros=total_zeros+size(Z{j},1);
end

y=zeros(nvar*nvar+total_zeros,1);

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

% p4 must be of norm one
p4 = x(end);
y(k)= norm(p4)^2 - 1.0;

k=k+1;


f    = FA0Aplus_f(f_h_P3_p4_inv(x,info),info.npredetermined,info.nvar,info.nlag,hirfs);

ib=k;
for j=1:n
    ie=ib+size(Z{j},1);  
    y(ib:ie-1)=Z{j}*f(:,j);
    ib=ie;
end
