function [St_old] = restriction_cmy_3shocks(Bdraw, Sigmadraw, Qdraw,fo_inv,fo_str2irfs,info)
% checking sign restrictions


    %BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];

    %structpara   = fo_inv(BSigmaQnew);

    %irfpara      = fo_str2irfs(structpara);


    hSigmadraw = chol(Sigmadraw);

    L0 = hSigmadraw'*Qdraw;%reshape(irfpara(1:info.nvar*info.nvar,:),info.nvar,info.nvar);





    St_old = all(info.Ss{1,1}*L0(:,1)>0)*all(info.Ss{2,1}*L0(:,2)>0)*all(info.Ss{3,1}*L0(:,3)>0);








