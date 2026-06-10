function [St_old] = restriction_baseline_sys(Bdraw, Sigmadraw, Qdraw,fo_BSigmaQ2OBJ,info)
% checking sign restrictions


BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];

objparanew = fo_BSigmaQ2OBJ(BSigmaQnew);

D = objparanew(1:info.nvar);

L0 = reshape(objparanew(info.nvar+1:info.nvar+info.nvar*info.nvar),info.nvar,info.nvar);


A0 = inv(L0)';


St_old =  (L0(1,1)>0)*(A0(1,1)>0)*(A0(2,1)<0)*(A0(2,1)>-4)*(A0(3,1)<0)*(A0(3,1)>-4)*(D(1,1)>0)*(L0(2,1)<0); % first shock and systematic






