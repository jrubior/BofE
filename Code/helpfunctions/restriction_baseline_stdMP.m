function [St_old] = restriction_baseline_stdMP(Bdraw, Sigmadraw, Qdraw,fo_BSigmaQ2OBJ,info)
% checking sign restrictions


BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];

objparanew = fo_BSigmaQ2OBJ(BSigmaQnew);

L0 = reshape(objparanew(info.nvar+1:info.nvar+info.nvar*info.nvar),info.nvar,info.nvar);


A0 = inv(L0)';



St_old = (A0(1,1)>0); % first shock and systematic






