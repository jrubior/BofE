function [Q,R,X] = DrawQ(n)
%
%
%

X=randn(n,n);
[Q,R]=qr(X);

for i=1:n
    if R(i,i) <0
        Q(:,i)=-Q(:,i);
        R(i,:)=-R(i,:);
    end
end
