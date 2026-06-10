function y = IRFToModulus2(x,info)
%
%  model:
%
%    y(t)'*A(0) = z(t)'*C + y(t-1)'*A(1) + ... + y(t-p)'*A(p) + epsilon(t)
%
%    y(t) - n x 1 endogenous variables
%    epsilon(t) - n x 1 exogenous shocks
%    z(t) - k x 1 exogenous or deterministic variables
%    A(i) - n x n
%    C - k x n   
%
%
%  IRF parameterization - L(0), L(1), ... L(p), C
%
%    A(0) = inv(L(0))
%    A(i) = inv(L(0))*(L(i)*A(0) - L(i-1)*A(1) - ... - L(1)*A(i-1))  0 < i <= p
%
%    x = vec([vec(L(0)'); ... vec(L(p)'); vec(C)])
%



n  = info.nvar;
p  = info.nlag;

n2 = n*n;
k  = numel(x)/n - n*(p+1);
L  = cell(p,1);
L0 = reshape(x(1:n2),n,n)';

OomegaTilde = diag(diag(L0))*diag(diag(L0));


A = cell([p,1]);
B = cell([p,1]);
for i=1:p
    L{i}=reshape(x(i*n2+1:(i+1)*n2),n,n)';
    X = L0\L{i};
    for j=1:i-1
        X = X - L{i-j}*A{j};
    end
    A{i} = L0\X;
    B{i} = A{i}*L0';
end

LTilde      = cell([p,1]);
L0Tilde     = L0/chol(OomegaTilde);
for i=1:p
    if i==1
        LTilde{i}=B{i}'*L0Tilde;
    else

        tmp = zeros(info.nvar);
        for kk=1:p
            if i-kk>0
                tmp = tmp + B{i}'*LTilde{i-kk};
            else
                tmp = tmp + B{i}'*L0Tilde;
            end

        end
        LTilde{i} = tmp;

    end

end


y=zeros(size(x,1)+info.nvar,1);

y(1:n2)=reshape(L0Tilde,n2,1);
for i=1:p
    y(i*n2+1:(i+1)*n2)=reshape(LTilde{i},n2,1);
end
c  =reshape(x(n2*(p+1)+info.nex:end),k,n);
y(n2*(p+1)+1:n2*(p+1)+info.nvar) = c/chol(OomegaTilde);
y(n2*(p+1)+info.nvar+1:end) = diag(OomegaTilde);
