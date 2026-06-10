function [b,M] = compactForecastNotation(A0,BETA,exog,tstart,hmax,y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

B = BETA(size(exog,2)+1:end,:);
c = BETA(1:size(exog,2),:);

n = size(A0,1);
L = size(B,1)/n;

B = permute(reshape(B',n,n,L),[2,1,3]);



% Create K matrices

K_h = nan(n,n,hmax);

K_h(:,:,1) = eye(n);

for i = 1:hmax-1
   sumterm = zeros(n,n);
   for j = 1:i
       if j <= L
           sumterm = sumterm + K_h(:,:,i+1-j)*B(:,:,j);
       end
   end
   K_h(:,:,i+1) = eye(n) + sumterm;
end


% Create N matrices

N_h_l = nan(n,n,hmax,L);

for i = 1
    for l = 1:L
        N_h_l(:,:,1,l) = B(:,:,l);
    end
end

for i = 2:hmax
   for l = 1:L
       sumterm = zeros(n,n);
       for j = 1:i-1
           if j <= L
           sumterm = sumterm + N_h_l(:,:,i-j,l)*B(:,:,j);
           end
       end
       if i+l-1 > L
           N_h_l(:,:,i,l) = sumterm;
       else
           N_h_l(:,:,i,l) = B(:,:,i+l-1) + sumterm;
       end
   end
end


% Create b matrix
b = [];

for h = 1:hmax
    
    sumterm = zeros(1,n);
    
    for l = 1:L
       sumterm  =  sumterm + y(tstart+1-l,:)*N_h_l(:,:,h,l);
    end
    
    b_tplush = c*K_h(:,:,h-1+1) + sumterm;
        
    b = [b b_tplush];

end

b = b';

% Create M matrix

M_i = nan(n,n,hmax);

M_i(:,:,1) = inv2(A0);
M = zeros(n*hmax);

for i = 1:hmax-1
   sumterm = zeros(n,n);
   for j = 1:i
       if j <= L
           sumterm = sumterm + M_i(:,:,i+1-j)*B(:,:,j);
       end
   end
   M_i(:,:,i+1) = sumterm;
   
end


Mi_reshaped = reshape(M_i,n,n*hmax);

for i = 1:hmax
 M(n*(i-1)+1:n*(i-1)+n,n*(i-1)+1:end) = Mi_reshaped(:,1:n*hmax-n*(i-1));
end


end

