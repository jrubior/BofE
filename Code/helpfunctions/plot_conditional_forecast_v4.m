function [fig_handles] = plot_conditional_forecast_v4( ...
    yf_save, yf_unc_save, yt, T, hmax, nvar, yCondition, VarNames, ...
    y_realized_future, dates_realized_future, dates_historical, dates_forecast, ...
    observableCondicionado_idx, show_realized, num_hist_periods, ...
    dates_fcast_xls, F_med_xls, F_p16_xls, F_p84_xls, ...
    tasas, vars_to_transform)
% PLOT_CONDITIONAL_FORECAST_V4
% - Versión con transformación a tasas (diferencias logarítmicas)
% - Mantiene la leyenda maestra exacta de la v3

% ==========================================================
% 1. Normalización de Fechas e Inputs
% ==========================================================
if istable(dates_historical), dates_historical = dates_historical{:,1}; end
if istable(dates_forecast),   dates_forecast   = dates_forecast{:,1};   end
if exist('dates_realized_future','var') && istable(dates_realized_future)
    dates_realized_future = dates_realized_future{:,1};
end
if exist('dates_fcast_xls','var') && istable(dates_fcast_xls)
    dates_fcast_xls = dates_fcast_xls{:,1};
end

if ~isdatetime(dates_historical), dates_historical = datetime(dates_historical); end
if ~isdatetime(dates_forecast),   dates_forecast   = datetime(dates_forecast);   end
if exist('dates_realized_future','var') && ~isempty(dates_realized_future) && ~isdatetime(dates_realized_future)
    dates_realized_future = datetime(dates_realized_future);
end
if exist('dates_fcast_xls','var') && ~isempty(dates_fcast_xls) && ~isdatetime(dates_fcast_xls)
    dates_fcast_xls = datetime(dates_fcast_xls);
end

use_xls = (nargin >= 19) && ~isempty(dates_fcast_xls) && ~isempty(F_med_xls) ...
          && ~isempty(F_p16_xls) && ~isempty(F_p84_xls);

% ==========================================================
% 2. Normalización de Datos a 3D [Time x Var x Draws]
% ==========================================================
if iscell(yf_save), yf_save_3D = cat(3, yf_save{:}); else, yf_save_3D = yf_save; end
if iscell(yf_unc_save), yf_unc_3D = cat(3, yf_unc_save{:}); else, yf_unc_3D = yf_unc_save; end

yf_save_3D = yf_save_3D(1:T+hmax, 1:nvar, :);
yf_unc_3D  = yf_unc_3D(1:T+hmax, 1:nvar, :);

% ==========================================================
% 3. Lógica de Transformación a Tasas (Diferencias)
% ==========================================================
transform_diff = @(data) [nan(1, size(data,2), size(data,3)); data(2:end,:,:) - data(1:end-1,:,:)];

yt_plot      = yt;
yf_cond_plot = yf_save_3D;
yf_unc_plot  = yf_unc_3D;
yCond_plot   = yCondition;
yReal_plot   = y_realized_future;

if tasas && ~isempty(vars_to_transform)
    for v = vars_to_transform
        yt_plot(2:end, v) = yt(2:end, v) - yt(1:end-1, v);
        yt_plot(1, v) = nan; 
        
        yf_cond_plot(:, v, :) = transform_diff(yf_save_3D(:, v, :));
        yf_unc_plot(:, v, :)  = transform_diff(yf_unc_3D(:, v, :));
        
        if use_xls
            full_xls_med = [yt(T,v); F_med_xls(:,v)];
            full_xls_p16 = [yt(T,v); F_p16_xls(:,v)];
            full_xls_p84 = [yt(T,v); F_p84_xls(:,v)];
            F_med_xls(:,v) = full_xls_med(2:end) - full_xls_med(1:end-1);
            F_p16_xls(:,v) = full_xls_p16(2:end) - full_xls_p16(1:end-1);
            F_p84_xls(:,v) = full_xls_p84(2:end) - full_xls_p84(1:end-1);
        end
        
        if ~isempty(yCondition)
            full_path = [yt(T,v); yCondition(:, v)];
            yCond_plot(:, v) = full_path(2:end) - full_path(1:end-1);
        end
        
        if ~isempty(y_realized_future)
            full_real = [yt(T,v); y_realized_future(:, v)];
            yReal_plot(:, v) = full_real(2:end) - full_real(1:end-1);
        end
    end
end

% ==========================================================
% 4. Estadísticos y Ventanas
% ==========================================================
cond_median = median(yf_cond_plot, 3, 'omitnan');
cond_p16    = prctile(yf_cond_plot, 16, 3);
cond_p84    = prctile(yf_cond_plot, 84, 3);
unc_median  = median(yf_unc_plot,  3, 'omitnan');

hist_idx = max(2, T - num_hist_periods + 1) : T;
time_hist_win = dates_historical(hist_idx);
yt_win = yt_plot(hist_idx, :);

dates_joined_model = [dates_historical(T); dates_forecast(:)];  
if use_xls
    H_common = 1 + min(hmax, length(dates_fcast_xls));
    dates_common = dates_joined_model(1:H_common);
    ia_model = (1:H_common)';
    ia_xls   = (1:H_common)';
else
    dates_common = dates_joined_model;
    ia_model = (1:length(dates_joined_model))';
end

% ==========================================================
% 5. Graficación
% ==========================================================
fig_handles = [];
max_per_fig = 9;
ncols = 3; nrows = 3;

% Definición de colores para consistencia
color_cond_CI = [0.8 0.8 1];
color_cond_med = 'b';
color_unc_med = [0.5 0.5 0.5];
color_xls_med = [0.4 0.4 0.4]; % Gris oscuro

for i = 1:nvar
    if mod(i-1, max_per_fig) == 0
        h_fig = figure('Name', sprintf('Forecast Page %d', ceil(i / max_per_fig)));
        set(h_fig, 'Position', [50, 50, 1200, 850]);
        fig_handles = [fig_handles; h_fig];
    end
    
    subplot(nrows, ncols, mod(i-1, max_per_fig) + 1);
    hold on;
    
    c_med = cond_median(T:end, i); c_med = c_med(ia_model);
    c_p16 = cond_p16(T:end, i);    c_p16 = c_p16(ia_model);
    c_p84 = cond_p84(T:end, i);    c_p84 = c_p84(ia_model);
    u_med = unc_median(T:end, i);   u_med = u_med(ia_model);
    
    % 1. Banda Condicional
    fill([dates_common; flipud(dates_common)], [c_p16; flipud(c_p84)], ...
         color_cond_CI, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    % 2. Mediana Condicional
    plot(dates_common, c_med, 'Color', color_cond_med, 'LineStyle', '-', 'LineWidth', 1.5);
    % 3. Historia
    plot(time_hist_win, yt_win(:, i), 'k-', 'LineWidth', 1.5);
    % 4. Incondicional
    plot(dates_common, u_med, ':', 'Color', color_unc_med, 'LineWidth', 1.5);
    
    % 5. Excel Forecast (Línea discontinua)
    if use_xls
        xls_med_j = [yt_plot(T,i); F_med_xls(:,i)];
        plot(dates_common, xls_med_j(ia_xls), '--', 'Color', color_xls_med, 'LineWidth', 1.5);
    end
    
    % 6. Condition Path (Verde)
    if ~isempty(observableCondicionado_idx) && i == observableCondicionado_idx
        cp_j = [yt_plot(T,i); yCond_plot(:, i)];
        plot(dates_common, cp_j(ia_model), 'g--', 'LineWidth', 1.8);
    end
    
    % 7. Realized (Rojo)
    if show_realized && ~isempty(yReal_plot)
        r_j = [yt_plot(T,i); yReal_plot(:, i)];
        plot(dates_common, r_j(ia_model), 'r-.', 'LineWidth', 1.5);
    end
    
    if tasas && ismember(i, vars_to_transform), yline(0, 'k:', 'Alpha', 0.2); end
    
    title(VarNames{i}); grid on;
    if tasas && ismember(i, vars_to_transform), ylabel('% Var Trim.'); else, ylabel('Valor'); end
    xlim([time_hist_win(1) dates_common(end)]);
    xtickangle(45);
end

% ==========================================================
% 6. LEYENDA MAESTRA SINCRONIZADA
% ==========================================================
try
    % Textos de la leyenda
    master_labels = {'Cond. 68% CI', 'History', 'Cond. Median', 'Uncond. Median'};
    if use_xls, master_labels = [master_labels, {'ForecastDraws Median'}]; end
    master_labels = [master_labels, {'Condition Path'}];
    if show_realized, master_labels = [master_labels, {'Realized'}]; end
    
    % Eje invisible para leyenda
    hL = axes('Visible','off', 'Position', [0.15 0.01 0.7 0.05], 'Parent', fig_handles(end));
    hold(hL, 'on');
    
    % Creamos los DUMMY PLOTS con los mismos estilos que arriba
    d = [];
    d(1) = fill(hL, NaN, NaN, color_cond_CI, 'EdgeColor', 'none'); % Banda azul
    d(2) = plot(hL, NaN, NaN, 'k-', 'LineWidth', 1.5);             % Historia negra
    d(3) = plot(hL, NaN, NaN, 'Color', color_cond_med, 'LineWidth', 1.5); % Mediana azul
    d(4) = plot(hL, NaN, NaN, ':', 'Color', color_unc_med, 'LineWidth', 1.5); % Incond. punteada
    
    idx = 5;
    if use_xls
        d(idx) = plot(hL, NaN, NaN, '--', 'Color', color_xls_med, 'LineWidth', 1.5); % Excel discontinua
        idx = idx + 1;
    end
    
    d(idx) = plot(hL, NaN, NaN, 'g--', 'LineWidth', 1.8); % Path verde
    idx = idx + 1;
    
    if show_realized
        d(idx) = plot(hL, NaN, NaN, 'r-.', 'LineWidth', 1.5); % Realized rojo
    end
    
    legend(hL, d, master_labels, 'Orientation','horizontal', 'Location','south', 'Box','off', 'FontSize', 9);
    hold(hL, 'off');
catch ME
    warning('Error en leyenda: %s', ME.message);
end

end