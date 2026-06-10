function y = SigmaQToOBJ(x,info)


nvar=info.nvar;
m=info.m;
nvar2=info.nvar^2;


Sigma=reshape(x(1:nvar2),nvar,nvar);
Q=reshape(x(nvar2+1:end),nvar,nvar);

D = diag(diag(chol(Sigma)'*Q));

L0= chol(Sigma)'*(Q/D);


y = nan(nvar2,1);


y(1:nvar)=reshape(diag(D),nvar,1);

vecL0 = L0(:);

index_diag = linspace(1,nvar2,nvar);

vecL0(index_diag) = [];

y(nvar+1:nvar2)=vecL0;






