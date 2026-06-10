function [St_old] = restriction_none(Bdraw, Sigmadraw, Qdraw,fo_inv,fo_str2irfs,info)
% checking sign restrictions


% BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];
% 
% objparanew = fo_BSigmaQ2OBJ(BSigmaQnew);
% 
% L0 = reshape(objparanew(info.nvar+1:info.nvar+info.nvar*info.nvar),info.nvar,info.nvar);
% 
% 
% % 
% % 
% % if L0(1,1)<0
% % 
% %     Qdraw(:,1)=-Qdraw(:,1);
% %     BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];
% % 
% %     objparanew = fo_BSigmaQ2OBJ(BSigmaQnew);
% % 
% %     L0 = reshape(objparanew(info.nvar+1:info.nvar+info.nvar*info.nvar),info.nvar,info.nvar);
% % 
% % end


St_old = 1;%L0(1,1)>0; % first shock and systematic






