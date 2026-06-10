function hFigs = plot_irf_h0_histograms(spec, variableNames, shocksToPlot, varargin)

p = inputParser;
addParameter(p, 'MatFolder','matfiles');
addParameter(p, 'CumulativeVarIdx', []);
addParameter(p, 'Scale', 1);
addParameter(p, 'ShockNames', {});
addParameter(p, 'VarOrderIdx', []);
addParameter(p, 'NumBins', 40);
addParameter(p, 'SavePNG', true);
addParameter(p, 'PngFolder', 'IRFs plots');

parse(p, varargin{:});
opt = p.Results;

matPath = fullfile(opt.MatFolder, [char(spec) '.mat']);
assert(exist(matPath,'file')==2, 'File not found: %s', matPath);
M = matfile(matPath);

wL    = whos(M, 'L');
wcumL = whos(M, 'cumL');
assert(~isempty(wL) && ~isempty(wcumL), 'MAT must contain variables L and cumL.');

dims = wL.size;
nVar   = dims(2);
nShock = dims(3);
nDraws = dims(4);

if isempty(opt.VarOrderIdx)
    ord = 1:nVar;
else
    ord = opt.VarOrderIdx(:)';
end

nVarPlot = numel(ord);
assert(numel(variableNames)==nVarPlot, ...
    'variableNames must have length %d.', nVarPlot);

shocksToPlot = shocksToPlot(:)';
assert(all(shocksToPlot>=1 & shocksToPlot<=nShock), 'Invalid shocksToPlot.');

cumIdxMat = unique(opt.CumulativeVarIdx(:)');
scale = opt.Scale;

if isempty(opt.ShockNames)
    shockTitles = arrayfun(@(k)sprintf('Shock %d',k), 1:nShock, 'UniformOutput', false);
else
    shockTitles = cellstr(opt.ShockNames);
end

gridRows = 4; %2; 
gridCols = 4; %3; 
perPage = gridRows * gridCols;
pagesPerShock = max(1, ceil(nVarPlot / perPage));

hFigs = gobjects(numel(shocksToPlot)*pagesPerShock,1);
hIdx = 0;

for si = 1:numel(shocksToPlot)
    k = shocksToPlot(si);

    for pgi = 1:pagesPerShock
        startVar = (pgi-1)*perPage + 1;
        endVar   = min(nVarPlot, pgi*perPage);
        nThis    = endVar - startVar + 1;

        hFig = figure('Name', sprintf('Histograms h=0 | %s | %s | page %d/%d', ...
            char(spec), shockTitles{k}, pgi, pagesPerShock));

        tl = tiledlayout(hFig, gridRows, gridCols, ...
            'Padding','compact', 'TileSpacing','compact');

        for j = 1:nThis
            jj = startVar + j - 1;
            v  = ord(jj);

            ax = nexttile(tl);
            hold(ax, 'on');

            if ismember(v, cumIdxMat)
                draws = squeeze(M.cumL(1, v, k, :));
            else
                draws = squeeze(M.L(1, v, k, :));
            end

            draws = draws(:) * scale;

            histogram(ax, draws, opt.NumBins, ...
                'Normalization', 'pdf');

            xline(ax, quantile(draws, 0.16), '--', 'LineWidth', 1.2);
            xline(ax, quantile(draws, 0.50), '-',  'LineWidth', 1.4);
            xline(ax, quantile(draws, 0.84), '--', 'LineWidth', 1.2);
            xline(ax, 0, 'k-.', 'LineWidth', 0.8);

            hold(ax, 'off');

            grid(ax, 'on');
            box(ax, 'off');
            title(ax, variableNames{jj}, 'FontWeight','normal');
            xlabel(ax, 'IRF at h=0');
            ylabel(ax, 'Density');
        end

        title(tl, sprintf('Distribution of IRFs at h=0: %s', shockTitles{k}), ...
            'FontWeight','normal');

        if opt.SavePNG
            specName = char(spec);
            pngDir = fullfile(opt.PngFolder, specName);
            ensure_dir(pngDir);

            set(hFig, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);

            suffix = sprintf('_h0_hist_shock_%d_page_%dof%d', ...
                k, pgi, pagesPerShock);

            print(hFig, fullfile(pngDir, [specName suffix '.png']), ...
                '-dpng', '-r150');
        end

        hIdx = hIdx + 1;
        hFigs(hIdx) = hFig;
    end
end

if hIdx < numel(hFigs)
    hFigs = hFigs(1:hIdx);
end

drawnow;
end