function y = IRFToModulus(x,info)



n  = info.nvar;
p  = info.nlag;

n2 = n*n;
k  = numel(x)/n - n*(p+1);

L0 = reshape(x(1:n2),n,n);

OomegaTilde = diag(diag(L0))*diag(diag(L0));


Ltilde      = cell([p,1]);
L0tilde     = L0/chol(OomegaTilde);
for i=1:p
    Ltilde{i}=reshape(x(i*n2+1:(i+1)*n2),n,n)/chol(OomegaTilde);
end

y=zeros(size(x,1)+info.nvar,1);

y(1:n2)=reshape(L0tilde,n2,1);
for i=1:p
    y(i*n2+1:(i+1)*n2)=reshape(Ltilde{i},n2,1);
end
c  =reshape(x(n2*(p+1)+info.nex:end),k,n);
y(n2*(p+1)+1:n2*(p+1)+info.nvar) = c/chol(OomegaTilde);
y(n2*(p+1)+info.nvar+1:end) = diag(OomegaTilde);





