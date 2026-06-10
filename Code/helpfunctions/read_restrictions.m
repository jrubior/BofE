function info = read_restrictions(info, excelPath, IRF_sheetSigns, IRF_sheetRanks, A0_sheetSigns, A0_sheetBounds)

    %% ------------------------------------
    % 1) Leer hoja de IRF SIGN restrictions
    % -------------------------------------

    C = readcell(excelPath, 'Sheet', IRF_sheetSigns);

    % Suponemos: fila 1 = nombres shocks, col 1 = nombres variables
    vals = C(2:end, 2:end);   % solo la parte numérica

    [nvar, nshocks] = size(vals);
    M = zeros(nvar, nshocks);

    % Convertir a numérico (texto → número, vacíos → 0)
    for i = 1:nvar
        for j = 1:nshocks
            v = vals{i,j};

            if isempty(v)
                M(i,j) = 0;

            elseif isnumeric(v)
                M(i,j) = v;

            elseif ischar(v) || isstring(v)
                vnum = str2double(strtrim(v));
                if isnan(vnum)
                    M(i,j) = 0;
                else
                    M(i,j) = vnum;
                end

            else
                M(i,j) = 0;
            end
        end
    end

    M(isnan(M)) = 0;   % por seguridad

    %% -------------------------------
    % 2) Construir info.Ss con signos
    % -------------------------------

    info.nshocks = nshocks;          
    info.Ss      = cell(info.nshocks,1);

    for sh = 1:nshocks
        col = M(:, sh);
        idx = find(col ~= 0);        % filas con restricción

        nr  = numel(idx);
        Ssh = zeros(nr, nvar);

        for r = 1:nr
            iVar = idx(r);
            Ssh(r, iVar) = col(iVar);
        end

        info.Ss{sh,1} = Ssh;
    end

    %% -------------------------------------
    % 3) Leer hoja de IRF RANK restrictions
    % --------------------------------------
    if nargin >= 4 && ~isempty(IRF_sheetRanks)
        
        Tr = readtable(excelPath, 'Sheet', IRF_sheetRanks);
        % columnas esperadas: shock, fila, var_1, valor_1, var_2, valor_2

        for sh = 1:info.nshocks

            Tr_sh = Tr(Tr.shock == sh, :);
            if isempty(Tr_sh)
                continue;
            end

            Ssh = info.Ss{sh,1};

            for k = 1:height(Tr_sh)
                newrow = zeros(1, nvar);

                v1   = Tr_sh.var_1(k);
                val1 = Tr_sh.valor_1(k);
                v2   = Tr_sh.var_2(k);
                val2 = Tr_sh.valor_2(k);

                if ~isnan(v1) && ~isnan(val1)
                    newrow(v1) = val1;
                end
                if ~isnan(v2) && ~isnan(val2)
                    newrow(v2) = val2;
                end

                Ssh = [Ssh; newrow];
            end

            info.Ss{sh,1} = Ssh;
        end
        
    end

    %% ------------------------------------
    % 4) Leer hoja de IRF ZERO restrictions
    % -------------------------------------
    
    info.ZeroRest = cell(info.nshocks,1);  % vectores auxiliares

    C = readcell(excelPath, 'Sheet', IRF_sheetSigns);
    vals = C(2:end, 2:end);  % Suponemos: fila 1 = nombres shocks, col 1 = nombres variables  
    zeroIdx = zeros(nvar, nshocks);
    
    for i = 1:nvar
        for j = 1:nshocks
            x = vals{i,j};
            if isnumeric(x) && isscalar(x) && x == 0
                zeroIdx(i,j) = 1; % posición de los ceros
            end
        end
    end


    for sh = 1:info.nshocks
        info.ZeroRest{sh,1} = zeros(size(info.Ss{sh,1},1),1);
        if sum(zeroIdx(:,sh))==0
            continue;
        end

        idx = find(zeroIdx(:,sh));
        ns1 = sum(zeroIdx(:,sh))*2; % (x2 por valor absoluto)
        newrows = zeros(ns1,nvar); 
        rows = 1:2;
        for j = 1:numel(idx)
            newrows(rows(1) , idx(j)) = -1;  
            newrows(rows(2) , idx(j)) =  1;  
            rows = rows+2;
            info.ZeroRest{sh,1} = [info.ZeroRest{sh,1}; info.IRFepsilon; info.IRFepsilon];
        end
        info.Ss{sh,1} = [ info.Ss{sh,1}; newrows];
    end    
     
    %% ------------------------------------
    % 5) Leer hoja de A0 SIGN restrictions
    % -------------------------------------
    % Las restricciones en A0 van asociadas a los shocks. Por ejemplo, las
    % restricciones asociadas al shock 'i' se deben imponer en la columna 
    % 'i' de A0. 
    if ~isempty(A0_sheetSigns)

        C = readcell(excelPath, 'Sheet', A0_sheetSigns);
        B = readcell(excelPath, 'Sheet', A0_sheetBounds);
        B = B(2:end, 2:end); 
        vals = C(2:end, 2:end);  % Suponemos: fila 1 = nombres shocks, col 1 = nombres variables
        [nvar, nshocks] = size(vals);
        M = zeros(nvar, nshocks);
        for i = 1:size(M,1)
            for j = 1:size(M,2)
                if strcmpi(string(C{i+1,j+1}),'N')
                    M(i,j) = 1;
                end
            end
        end
        info.A0e = M; % Cada columna representa un shock. La fila del 1 indica el coeficiente de la variable por la que se va a normalizar.
        
        info.A0restrictionsSIGN = sum(any(M ~= 0, 1)); % número de shocks en los cuales restringimos coeficientes estructurales
        info.A0Ss = cell(info.A0restrictionsSIGN,1);
        info.A0Bounds = cell(info.A0restrictionsSIGN,1);

        for i = 1:info.A0restrictionsSIGN
            % Columna i de A0
            col = vals(:,i);
            ns1 = sum(cellfun(@(x) isnumeric(x) && isscalar(x) && abs(x)==1, col)) + 1;
            info.A0Ss{i} = zeros(ns1,nvar);
            info.A0Bounds{i} = zeros(ns1,1); % Desigualdades 
            
            idx = find(cellfun(@(x) (isnumeric(x) && isscalar(x) && abs(x)==1) || ((ischar(x) || isstring(x)) && strcmp(string(x),'N')), col));
            val = cellfun(@(x) double(strcmp(string(x),'N')) + double(isnumeric(x) && isscalar(x))*x, col(idx));
            
            idx_b = find(~cellfun(@ismissing,B(:,i)));

            % Signos A0
            for j =  1:ns1
                info.A0Ss{i}(j,idx(j)) = val(j);   
            end
            
            % Bounds (if any)
            if ~isempty(idx_b)
                val_b = B(idx_b,i);
                [~, idx_b] = ismember(idx_b,idx);
                for j = 1:numel(idx_b)
                    info.A0Bounds{i}(idx_b(j),1) = val_b{j}; % coeficiente > val
                end
            end
        end

    end

    %% ------------------------------------------
    % 6) Leer hoja de A0 QUASI-ZERO restrictions
    % -------------------------------------------
    if ~isempty(A0_sheetSigns)
        
        C = readcell(excelPath, 'Sheet', A0_sheetSigns);
        vals = C(2:end, 2:end);  % Suponemos: fila 1 = nombres shocks, col 1 = nombres variables
        [nvar, nshocks] = size(vals);
        zeroIdx = zeros(nvar, nshocks);  
        for i = 1:nvar
            for j = 1:nshocks
                x = vals{i,j};
                if isnumeric(x) && isscalar(x) && x == 0
                    zeroIdx(i,j) = 1; % posición de los ceros
                end
            end
        end

        info.A0restrictionsZERO = info.A0restrictionsSIGN;  
        info.A0Zs = cell(info.A0restrictionsZERO,1);
        info.A0eZEROS = zeros(nvar,info.A0restrictionsZERO); % columnas indican shocks, filas indican la variable por la que normalizamos
        shocks = find(any(cellfun(@(x) isnumeric(x) && isscalar(x) && x==0, vals),1));
        info.A0eZEROS(:,shocks) = info.A0e(:,shocks);

        for i = 1:numel(shocks)
            idx = find(zeroIdx(:,shocks(i)));
            ns1 = sum(zeroIdx(:,shocks(i)))*2; % (x2 por valor absoluto)
            info.A0Zs{shocks(i)} = zeros(ns1,nvar); 
            rows = 1:2;
            for j = 1:numel(idx)
                info.A0Zs{shocks(i)}(rows(1) , idx(j)) = -1;  
                info.A0Zs{shocks(i)}(rows(2) , idx(j)) =  1;  
                rows = rows+2;
            end
        end

    end

end

