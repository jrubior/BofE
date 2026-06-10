function [St_old] = restriction_baseline_uhlig(Bdraw, Sigmadraw, Qdraw,fo_BSigmaQ2OBJ,info)
% checking sign restrictions


BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];

objparanew = fo_BSigmaQ2OBJ(BSigmaQnew);

L0 = reshape(objparanew(info.nvar+1:info.nvar+info.nvar*info.nvar),info.nvar,info.nvar);




if (L0(1,1)<0)*(L0(3,1)>0)==1

    Qdraw(:,1)=-Qdraw(:,1);
    BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];

    objparanew = fo_BSigmaQ2OBJ(BSigmaQnew);

    L0 = reshape(objparanew(info.nvar+1:info.nvar+info.nvar*info.nvar),info.nvar,info.nvar);

end


St_old =  (L0(1,1)>0)*(L0(3,1)<0); % first shock and systematic






