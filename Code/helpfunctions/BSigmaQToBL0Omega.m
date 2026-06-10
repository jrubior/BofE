function y = BSigmaQToBL0Omega(x,info)


nvar=info.nvar;
m=info.m;
nvar2=info.nvar^2;

B=reshape(x(1:m*nvar),m,nvar);
Sigma=reshape(x(m*nvar+1:(m+nvar)*nvar),nvar,nvar);
Q=reshape(x((m+nvar)*nvar+1:end),nvar,nvar);

Omega = diag(diag(chol(Sigma)'*Q))^2;

D = chol(Omega);

L0= chol(Sigma)'*(Q/D);

keyboard
if info.nex==1
y = nan(nvar2+nvar2*info.nlag+nvar,1);
else
y = nan(nvar2+nvar2*info.nlag,1);
end

%y(1:nvar)=reshape(diag(D),nvar,1);
y(1:nvar)=reshape(diag(Omega),nvar,1);

vecL0 = L0(:);

index_diag = linspace(1,nvar2,nvar);

vecL0(index_diag) = [];

y(nvar+1:nvar2)=vecL0;


Lplus = nan(m,nvar);

for ell =1:info.nlag
 

        tmp = zeros(info.nvar);

        for k=1:min(ell,info.nlag)
            if ell-k==0
                tmp = tmp + B(info.nex+1+(ell-1)*info.nvar:info.nex+ell*info.nvar,:)'*L0;
            else
                tmp =  tmp+ B(info.nex+1+(k-1)*info.nvar:info.nex+k*info.nvar,:)'*Lplus(1+((ell-k)-1)*info.nvar:(ell-k)*info.nvar,:);
            end
        end
        Lplus(1+(ell-1)*info.nvar:ell*info.nvar,:)=tmp;
end



for ell=1:info.nlag
    y(ell*nvar2+1:(ell+1)*nvar2)=reshape(Lplus(1+(ell-1)*info.nvar:ell*info.nvar,:),nvar2,1);
end

if info.nex==1
c = inv(D)'*Q'*inv(chol(Sigma))'*(B(1,:)');
y((info.nlag+1)*nvar2+1:end)=c;

end

