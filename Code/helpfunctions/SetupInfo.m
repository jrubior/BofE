function info = SetupInfo(nvar,npredetermined,nlag, horizons, h)


NS        = 1 + numel(horizons);

%==========================================================================
%% identification Ss Zs
%==========================================================================


% NO sign restrictions
info.S = cell(nvar,1);
for ii=1:nvar
    info.S{ii}=zeros(0,nvar*NS);
end


% NO zero restrictions
info.Z=cell(nvar,1);
for i=1:nvar
    info.Z{i}=zeros(0,nvar*NS);
end



% computes the dimension of space in which the spheres live that map into
% the orthogonal matrices satisfying the restriction.  
dim=0;
for i=1:nvar
    dim=dim+nvar-(i-1+size(info.Z{i},1));
end

dimunc=0;
for i=1:nvar
    dimunc=dimunc+nvar-(i-1);
end

% number zero restrictions
nzeros=0;
for i=1:nvar
    nzeros=nzeros+size(info.Z{i},1);
end


nzerosunc=0;
for i=1:nvar
    nzerosunc=nzerosunc+size(info.Z{i},1);
end

% gets random W
W=cell(nvar,1);
for j=1:nvar
    W{j}=randn(nvar-(j-1+size(info.Z{j},1)),nvar);
end

% gets random W
Wunc=cell(nvar,1);
for j=1:nvar
    Wunc{j}=randn(nvar-(j-1),nvar);
end



% info
info.nvar=nvar;
info.m=npredetermined;
info.h=h;
info.dim=dim;
info.nzeros=nzeros;
info.W=W;
info.dimunc=dimunc;
info.nzerosunc=nzerosunc;
info.Wunc=Wunc;
info.e = eye(nvar);

info.nlag     = nlag;
info.horizons = horizons;
info.ZF       = @(x,y)ZF(x,y); % Z times F(A_{0},A_{+})
info.nex      = npredetermined-nvar*nlag;
info.nlag     = nlag;
info.nvar     = nvar;