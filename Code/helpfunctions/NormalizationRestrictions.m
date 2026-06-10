function r = NormalizationRestrictions(x, info)
%
%
%
n=info.nvar;
n2=n*n;
vecL0= x(1+n:n2+n);
L0 = reshape(vecL0,n,n);

r=zeros(n,1); %number of restrictions


% diagonal of L0 must be of norm one
for i=1:n
    r(i,1)=L0(i,i)^2 - 1.0;
end