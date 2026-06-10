function y = ModulusToIRF(x,info)


n  = info.nvar;
p  = info.nlag;

n2 = n*n;
k  = (numel(x)-n)/n - n*(p+1);

L0tilde    = reshape(x(1:n2),n,n);

OomegaTilde = diag(x(end-n+1:end));


L     = cell([p,1]);
L0    = L0tilde*chol(OomegaTilde);
for i=1:p
    L{i}=reshape(x(i*n2+1:(i+1)*n2),n,n)*chol(OomegaTilde);
end

y=zeros(size(x,1)-info.nvar,1);

y(1:n2)=reshape(L0,n2,1);
for i=1:p
    y(i*n2+1:(i+1)*n2)=reshape(L{i},n2,1);
end

ctilde  =reshape(x(n2*(p+1)+info.nex:n2*(p+1)+info.nvar*info.nex),k,n);

y(n2*(p+1)+info.nex:end) = ctilde*chol(OomegaTilde);



