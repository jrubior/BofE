function [q,r] = qr_unique(x)

[q,r] = qr(x);

for i=1:size(x,1)

    if r(i,i)<0
        q(:,i)= -q(:,i);
    end

end

end