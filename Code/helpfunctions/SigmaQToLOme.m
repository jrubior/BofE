function y = SigmaQToLOme(x,info)


nvar=info.nvar;
m=info.m;
nvar2=info.nvar^2;


Sigma=reshape(x(1:nvar2),nvar,nvar);
Q=reshape(x(nvar2+1:end),nvar,nvar);

Ome = diag(diag(h_tilde(Sigma,info)'*Q))^2;% diag(diag(h_tilde(Sigmam,info)'*Q))*diag(diag(h_tilde(Sigma,info)'*Q));
D = h_tilde(Ome,info);
L0= h_tilde(Sigma,info)'*(Q/D);

y = nan(nvar2,1);


y(1:nvar)=reshape(diag(Ome),nvar,1);

vecL0 = L0(:);

% index_diag = linspace(1,nvar2,nvar);
% 
% vecL0(index_diag) = [];

y(nvar+1:nvar2+nvar)=vecL0;





