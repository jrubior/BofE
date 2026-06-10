function SRvec = SRvec(Bdraw, Sigmadraw, Qdraw, info)

    % checking sign restrictions
    nvar = info.nvar;
    npredetermined = info.m;
    hSigmadraw = chol(Sigmadraw);
    
    L0 = hSigmadraw'*Qdraw;

    % structural coeficients
    BSigmaQ = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];
    structpara = f_h_inv(BSigmaQ,info);
    n2=nvar*nvar;
    A0=reshape(structpara(1:n2),nvar,nvar);

    SRvec = [];
    count = 1;

    % --- Restricciones de SIGNO en coeficientes estructurales contemporáneos: A0
    if isfield(info,'A0restrictionsSIGN') && ~isempty(info.A0restrictionsSIGN)
        for i=1:info.A0restrictionsSIGN
            n_A0restrictions_i = size(info.A0Ss{i,1},1);
            for j = 1:n_A0restrictions_i
                coef_denominador = info.A0e(:,i)'*A0*info.e(:,i); 
                coef_numerador = info.A0Ss{i,1}(j,:)*A0*info.e(:,i);
                SRvec(count,1) = (coef_numerador/coef_denominador) + info.A0Bounds{i}(j,1);
                count = count + 1;
            end
        end
    end

    % --- Restricciones de CEROS en coeficientes estructurales contemporaneos: A0
    if  isfield(info,'A0restrictionsZERO') && ~isempty(info.A0restrictionsZERO)
        for i = 1:info.A0restrictionsZERO
            n_A0restrictions_i = size(info.A0Zs{i,1},1);
            if n_A0restrictions_i > 0
                for j = 1:n_A0restrictions_i  
                    coef_denominador = info.A0eZEROS(:,i)'*A0*info.e(:,i); 
                    coef_numerador = info.A0Zs{i,1}(j,:)*A0*info.e(:,i);             
                    SRvec(count,1) = (coef_numerador/coef_denominador) + info.A0epsilon;
                    count = count + 1;
                end
            end
        end
    end

    % --- Contemporáneo: L0
    for i=1:info.nshocks
        n_restrictions_i = size(info.Ss{i,1},1);
        for j=1:n_restrictions_i
            SRvec(count,1) = info.Ss{i,1}(j,:)*L0*info.e(:,i) + info.ZeroRest{i,1}(j); % modificación para quasi-zero restrictions:   |impacto| < epsilon => (- impacto + epsilon > 0) & (impacto + epsilon > 0)
            count = count + 1;
        end
    end

    % --- (2) Restricciones dinámicas (si existen)
    if isfield(info,'dyn') && ~isempty(info.dyn) && isfield(info.dyn,'restr') && ~isempty(info.dyn.restr)
    
        horizons = info.dyn.horizons;
        IRFstack = IRF_horizons(structpara, nvar, info.nlag, info.pred, info.nex, horizons);
        restr = info.dyn.restr;

        if ~isstruct(restr), error('info.dyn.restr debe ser struct'); end
        if numel(restr) == 1, restr = restr(:); end
    
        for rr = 1:numel(restr)
            v  = restr(rr).v;      % variable (posición, fila de L_k)
            sh = restr(rr).s;      % shock (posición, columna de L_k)
            sg = restr(rr).sign;   % +1 => signo positivo, -1 => signo negativo
            
            % Verificación de restricciones
            for k = 1:numel(horizons)
                rows = (k-1)*nvar + (1:nvar);
                Rh  = IRFstack(rows,:); 
                resp = Rh(v, sh);
                SRvec(count:count+size(sg,2)-1,1) = sg' .* resp; % imponer resp > 0 si sg=+1, resp <0 si sg=-1
                count = count + size(sg,2);
            end
        end
    end

end