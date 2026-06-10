function [St_old] = restriction_cmy(Bdraw, Sigmadraw, Qdraw,fo_inv,fo_str2irfs,info)
% checking sign restrictions


    %BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];

    %structpara   = fo_inv(BSigmaQnew);

    %irfpara      = fo_str2irfs(structpara);


    hSigmadraw = chol(Sigmadraw);

    L0 = hSigmadraw'*Qdraw;%reshape(irfpara(1:info.nvar*info.nvar,:),info.nvar,info.nvar);


St_old=1;

for j=1:info.nshocks
    St_old = St_old*all(info.Ss{j,1}*L0(:,j)>0); 
end







