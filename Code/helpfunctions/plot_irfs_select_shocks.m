function hFigs = plot_irfs_select_shocks(spec, variableNames, shocksToPlot, varargin)
% Plot IRFs in 4x4 grids, one figure per selected shock (auto-paginate if >16 vars).
%
% Expects MAT file at MatFolder/spec.mat with variables:
%   L      : [H x nVar x nShock x nDraws]  (non-cumulative IRFs)
%   cumL   : [H x nVar x nShock x nDraws]  (cumulative IRFs)
%   horizon: scalar with H == horizon+1
%
% New:
%   Saves IRFs to Excel, one sheet per shock.

% ----- Parse inputs -----
p = inputParser;
addParameter(p, 'MatFolder','matfiles');
addParameter(p, 'CumulativeVarIdx', []);
addParameter(p, 'YLimByVar', []);
addParameter(p, 'Scale', 1);
addParameter(p, 'SavePNG', true);
addParameter(p, 'SaveEPS', true);
addParameter(p, 'PngFolder', 'IRF plots');
addParameter(p, 'EpsFolder', 'epsfiles');
addParameter(p, 'ShockNames', {});
addParameter(p, 'VarOrderIdx', []);

% NUEVOS parámetros para Excel
addParameter(p, 'SaveExcel', true);
addParameter(p, 'ExcelFolder', 'IRF excels');
addParameter(p, 'ExcelFile', '');   % si está vacío, usa <spec>_IRFs.xlsx

parse(p, varargin{:});
opt = p.Results;

% ----- Load (MATFILE, no load of big arrays) -----
matPath = fullfile(opt.MatFolder, [char(spec) '.mat']);
assert(exist(matPath,'file')==2, 'File not found: %s', matPath);
M = matfile(matPath);

% Validate required variables exist
wL    = whos(M, 'L');
wcumL = whos(M, 'cumL');
whor  = whos(M, 'horizon');
assert(~isempty(wL) && ~isempty(wcumL) && ~isempty(whor), ...
    'MAT must contain variables L, cumL, horizon.');

dims = wL.size; % [H x nVar x nShock x nDraws]
H      = dims(1);
nVar   = dims(2);
nShock = dims(3);
nDraws = dims(4);

x = 0:M.horizon;
assert(H==numel(x), 'Mismatch: size(L,1)=%d but horizon+1=%d.', H, numel(x));

% ----- Optional: reorder variables according to VarOrderIdx -----
if isempty(opt.VarOrderIdx)
    ord = 1:nVar;  % keep original
else
    ord = opt.VarOrderIdx(:)';  % row vector
    assert(all(ord>=1 & ord<=nVar), 'VarOrderIdx has invalid indices.');
end
nVarPlot = numel(ord);

assert(numel(variableNames)==nVarPlot, ...
    'variableNames must have length %d (after VarOrderIdx).', nVarPlot);

% shocks
shocksToPlot = shocksToPlot(:)';
assert(all(shocksToPlot>=1 & shocksToPlot<=nShock), 'Invalid shocksToPlot indices.');

% ----- Plot settings -----
gridRows = 4; gridCols = 4; perPage = gridRows*gridCols;
% gridRows = 2; gridCols = 3; perPage = gridRows*gridCols;
pagesPerShock = max(1, ceil(nVarPlot / perPage));
fsAxis=10; fsLabel=10; fsTitle=11; lineW=1.4;
scale = opt.Scale;

% Cum vars & quantile probs
cumIdxMat = unique(opt.CumulativeVarIdx(:)');  % indices in MAT order
probs = [0.16 0.50 0.84];

% Shock names (fallbacks)
if isempty(opt.ShockNames)
    shockTitles = arrayfun(@(k)sprintf('Shock %d',k), 1:nShock, 'UniformOutput', false);
else
    shockTitles = cellstr(opt.ShockNames);
    if numel(shockTitles) < nShock
        shockTitles(end+1:nShock) = arrayfun(@(k)sprintf('Shock %d',k), ...
            (numel(shockTitles)+1):nShock, 'UniformOutput', false);
    end
end

% ----- Excel setup -----

if opt.SaveExcel
    specName = char(spec);
    excelDir = opt.ExcelFolder;
    ensure_dir(excelDir);

    if isempty(opt.ExcelFile)
        excelPath = fullfile(excelDir, [specName '_IRFs.xlsx']);
    else
        excelPath = fullfile(excelDir, opt.ExcelFile);
    end

    % Borra archivo previo para evitar hojas viejas si vuelves a correr
    if exist(excelPath, 'file') == 2
        delete(excelPath);
    end
end

% ----- Plot -----
hFigs = gobjects(numel(shocksToPlot)*pagesPerShock,1);
hIdx = 0;

for si = 1:numel(shocksToPlot)
    k = shocksToPlot(si);

    % ============================================================
    % 1) CALCULAR CUANTILES PARA TODAS LAS VARIABLES DEL SHOCK k
    %    (esto sirve tanto para graficar como para exportar Excel)
    % ============================================================
    vMatAll = ord;                    % orden final de ploteo
    labelsAll = variableNames(:)';    % etiquetas en orden final
    isCumAll = ismember(vMatAll, cumIdxMat);

    vMinAll = min(vMatAll);
    vMaxAll = max(vMatAll);
    vBlockAll = vMinAll:vMaxAll;
    nBlockAll = numel(vBlockAll);

    LblockAll    = squeeze(M.L(:,    vBlockAll, k, :));   % [H x nBlock x nDraws]
    cumLblockAll = squeeze(M.cumL(:, vBlockAll, k, :));   % [H x nBlock x nDraws]

    if nBlockAll == 1
        LblockAll    = reshape(LblockAll,    H, 1, nDraws);
        cumLblockAll = reshape(cumLblockAll, H, 1, nDraws);
    end

    posAll = vMatAll - vMinAll + 1;

    q16_shock = nan(H, nVarPlot);
    q50_shock = nan(H, nVarPlot);
    q84_shock = nan(H, nVarPlot);

    idxNon = find(~isCumAll);
    if ~isempty(idxNon)
        posNon = posAll(idxNon);
        dataNon = LblockAll(:, posNon, :);      % [H x nNon x nDraws]
        tmp = permute(dataNon, [3 1 2]);        % [nDraws x H x nNon]
        Q  = quantile(tmp, probs, 1);           % [3 x H x nNon]
        Q  = permute(Q, [2 3 1]);               % [H x nNon x 3]
        q16_shock(:, idxNon) = Q(:,:,1);
        q50_shock(:, idxNon) = Q(:,:,2);
        q84_shock(:, idxNon) = Q(:,:,3);
    end

    idxCum = find(isCumAll);
    if ~isempty(idxCum)
        posCum = posAll(idxCum);
        dataCum = cumLblockAll(:, posCum, :);   % [H x nCum x nDraws]
        tmp = permute(dataCum, [3 1 2]);        % [nDraws x H x nCum]
        Q  = quantile(tmp, probs, 1);           % [3 x H x nCum]
        Q  = permute(Q, [2 3 1]);               % [H x nCum x 3]
        q16_shock(:, idxCum) = Q(:,:,1);
        q50_shock(:, idxCum) = Q(:,:,2);
        q84_shock(:, idxCum) = Q(:,:,3);
    end

    % aplicar escala también a lo exportado
    q16_shock = q16_shock * scale;
    q50_shock = q50_shock * scale;
    q84_shock = q84_shock * scale;

    % ============================================================
    % 2) EXPORTAR A EXCEL: una hoja por shock
    % ============================================================
    if opt.SaveExcel
        excelData = table(x(:), 'VariableNames', {'Horizon'});

        for j = 1:nVarPlot
            baseName = matlab.lang.makeValidName(labelsAll{j});
            excelData.(sprintf('%s_q16', baseName)) = q16_shock(:, j);
            excelData.(sprintf('%s_q50', baseName)) = q50_shock(:, j);
            excelData.(sprintf('%s_q84', baseName)) = q84_shock(:, j);
        end

        sheetName = sanitize_sheet_name(shockTitles{k}, k);
        writetable(excelData, excelPath, 'Sheet', sheetName, 'WriteMode', 'overwritesheet');
    end

    % ============================================================
    % 3) GRAFICAR por páginas usando los cuantiles ya calculados
    % ============================================================
    for pgi = 1:pagesPerShock
        startVar = (pgi-1)*perPage + 1;
        endVar   = min(nVarPlot, pgi*perPage);
        nThis    = endVar - startVar + 1;

        hFig = figure('Name', sprintf('IRFs %s | %s | page %d/%d', ...
                      char(spec), shockTitles{k}, pgi, pagesPerShock));
        tl = tiledlayout(hFig, gridRows, gridCols, 'Padding','compact', 'TileSpacing','compact');

        labelsPage = labelsAll(startVar:endVar);

        for j = 1:nThis
            jj = startVar + j - 1;  % índice global dentro del shock

            ax = nexttile(tl);
            hold(ax, 'on');

            q16 = (q16_shock(:,jj))';
            q50 = (q50_shock(:,jj))';
            q84 = (q84_shock(:,jj))';

            plot(ax, x, q16, '--', 'LineWidth', lineW);
            plot(ax, x, q84, '--', 'LineWidth', lineW);
            plot(ax, x, q50, '-',  'LineWidth', lineW);
            plot(ax, x, zeros(size(x)), 'k-.', 'LineWidth', 0.5);

            hold(ax, 'off');

            grid(ax, 'on'); box(ax, 'off');
            set(ax, 'FontSize', fsAxis, 'LineWidth', 1.0);
            xlabel(ax, 'Quarters','FontSize',fsLabel);
            ylabel(ax, labelsPage{j},'FontSize',fsLabel);
            title(ax, labelsPage{j}, 'FontSize', fsTitle, 'FontWeight','normal');

            if ~isempty(opt.YLimByVar)
                yl = opt.YLimByVar(jj,:);   % índices LOCALES en orden de ploteo
                if numel(yl)==2 && all(isfinite(yl)), ylim(ax, yl); end
            end
        end

        if verLessThan('matlab','9.9')
            sgtitle(sprintf('IRFs: %s', shockTitles{k}), 'FontWeight','normal');
        else
            title(tl, sprintf('IRFs: %s', shockTitles{k}), 'FontWeight','normal');
        end

        if opt.SavePNG || opt.SaveEPS
            specName = char(spec);

            pngDir = fullfile(opt.PngFolder, specName);
            if opt.SavePNG
                ensure_dir(pngDir);
            end

            epsDir = fullfile(opt.EpsFolder, specName);
            if opt.SaveEPS
                ensure_dir(epsDir);
            end

            set(hFig, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);
            set(hFig, 'PaperUnits', 'inches');
            set(hFig, 'PaperPosition', [0 0 12 8]);
            set(hFig, 'PaperSize', [12 8]);
            set(hFig, 'PaperPositionMode', 'manual');

            suffix = sprintf('_shock_%d_page_%dof%d', k, pgi, pagesPerShock);

            if opt.SavePNG
                print(hFig, fullfile(pngDir, [specName suffix '.png']), '-dpng', '-r150');
            end
            if opt.SaveEPS
                print(hFig, fullfile(epsDir, [specName suffix '.eps']), '-depsc');
            end
        end

        hIdx = hIdx + 1;
        hFigs(hIdx) = hFig;
    end
end

if hIdx < numel(hFigs), hFigs = hFigs(1:hIdx); end
drawnow;
end


% =========================
% Helpers
% =========================
function sheetName = sanitize_sheet_name(name, k)
    sheetName = regexprep(char(name), '[:\\/?*\[\]]', '_'); % caracteres no válidos
    if isempty(strtrim(sheetName))
        sheetName = sprintf('Shock_%d', k);
    end
    if numel(sheetName) > 31
        sheetName = sheetName(1:31); % límite de Excel
    end
end

function ensure_dir(d)
    if exist(d,'dir')~=7, mkdir(d); end
end
