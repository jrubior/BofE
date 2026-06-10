function [fig_handles] = plot_conditional_forecast_v3( ...
    yf_save, yf_unc_save, yt, T, hmax, nvar, yCondition, VarNames, ...
    y_realized_future, dates_realized_future, dates_historical, dates_forecast, ...
    observableCondicionado_idx, show_realized, num_hist_periods, ...
    dates_fcast_xls, F_med_xls, F_p16_xls, F_p84_xls)

% PLOT_CONDITIONAL_FORECAST_V3
% - Condicional vs Incondicional (del run condicional) + Realized + Condition path
% - + Forecast_draws (incondicional corrida separada) resumido desde Excel (median/p16/p84)
% - Comparación SIEMPRE sobre el menor común por FECHAS (intersect), robusto a hmax ~= horizon_fcast.
% ==========================================================
% Normalizar inputs de fechas (por si llegan como table)
% ==========================================================
if istable(dates_historical), dates_historical = dates_historical{:,1}; end
if istable(dates_forecast),   dates_forecast   = dates_forecast{:,1};   end

if exist('dates_realized_future','var') && istable(dates_realized_future)
    dates_realized_future = dates_realized_future{:,1};
end
if exist('dates_fcast_xls','var') && istable(dates_fcast_xls)
    dates_fcast_xls = dates_fcast_xls{:,1};
end

% Convertir a datetime si vienen como strings/cellstr
if ~isdatetime(dates_historical), dates_historical = datetime(dates_historical); end
if ~isdatetime(dates_forecast),   dates_forecast   = datetime(dates_forecast);   end
if exist('dates_realized_future','var') && ~isempty(dates_realized_future) && ~isdatetime(dates_realized_future)
    dates_realized_future = datetime(dates_realized_future);
end
if exist('dates_fcast_xls','var') && ~isempty(dates_fcast_xls) && ~isdatetime(dates_fcast_xls)
    dates_fcast_xls = datetime(dates_fcast_xls);
end

%% --- Validaciones base ---
if length(dates_historical) ~= T
    error('La longitud de dates_historical (%d) no coincide con T (%d).', length(dates_historical), T);
end
if length(dates_forecast) ~= hmax
    error('La longitud de dates_forecast (%d) no coincide con hmax (%d).', length(dates_forecast), hmax);
end

use_xls = (nargin >= 19) && ~isempty(dates_fcast_xls) && ~isempty(F_med_xls) ...
          && ~isempty(F_p16_xls) && ~isempty(F_p84_xls);

if use_xls
    if size(F_med_xls,2) ~= nvar
        error('F_med_xls tiene %d variables, pero nvar=%d.', size(F_med_xls,2), nvar);
    end
    if size(F_p16_xls,2) ~= nvar || size(F_p84_xls,2) ~= nvar
        error('F_p16_xls/F_p84_xls no coinciden con nvar=%d.', nvar);
    end
    if size(F_med_xls,1) ~= length(dates_fcast_xls) || size(F_p16_xls,1) ~= length(dates_fcast_xls) || size(F_p84_xls,1) ~= length(dates_fcast_xls)
        error('Dimensiones Excel inconsistentes: filas de F_*_xls deben coincidir con length(dates_fcast_xls).');
    end
end

% ==========================================================
% Normalizar yf_save / yf_unc_save a formato 3D numérico
% Acepta: cell{ndraws} con (T+hmax x nvar) o array 3D directo
% ==========================================================
if iscell(yf_save)
    yf_save_3D = cat(3, yf_save{:});
else
    yf_save_3D = yf_save;
end

if iscell(yf_unc_save)
    yf_unc_3D = cat(3, yf_unc_save{:});
else
    yf_unc_3D = yf_unc_save;
end

% Validaciones básicas
if ~isnumeric(yf_save_3D) || ~isnumeric(yf_unc_3D)
    error('yf_save/yf_unc_save deben ser numéricos o celdas de numéricos.');
end
if ~isreal(yf_save_3D), yf_save_3D = real(yf_save_3D); end
if ~isreal(yf_unc_3D),  yf_unc_3D  = real(yf_unc_3D);  end

%% --- Ventana histórica a mostrar ---
hist_len = num_hist_periods;
hist_idx = max(1, T - hist_len + 1) : T;
time_hist_win = dates_historical(hist_idx);
yt_win = yt(hist_idx, :);

%% --- Cargar draws y calcular percentiles ---
% num_draws = length(yf_save);
% 
% all_forecasts_cond = nan(T + hmax, nvar, num_draws);
% all_forecasts_unc  = nan(T + hmax, nvar, num_draws);
% 
% for d = 1:num_draws
%     if ~isempty(yf_save{d}) && size(yf_save{d},1) == T + hmax && size(yf_save{d},2) == nvar
%         all_forecasts_cond(:,:,d) = yf_save{d};
%     end
%     if ~isempty(yf_unc_save{d}) && size(yf_unc_save{d},1) == T + hmax && size(yf_unc_save{d},2) == nvar
%         all_forecasts_unc(:,:,d) = yf_unc_save{d};
%     end
% end

% Ahora ya están en 3D: (T+hmax) x nvar x ndraws
all_forecasts_cond = yf_save_3D;
all_forecasts_unc  = yf_unc_3D;

% Si por alguna razón vienen con más variables, recorta a nvar:
all_forecasts_cond = all_forecasts_cond(:, 1:min(nvar,size(all_forecasts_cond,2)), :);
all_forecasts_unc  = all_forecasts_unc(:, 1:min(nvar,size(all_forecasts_unc,2)), :);

% Validar tamaño temporal mínimo
if size(all_forecasts_cond,1) < T + hmax || size(all_forecasts_unc,1) < T + hmax
    error('yf_save/yf_unc_save tienen menos filas que T+hmax. Cond=%d, Unc=%d, esperado=%d', ...
        size(all_forecasts_cond,1), size(all_forecasts_unc,1), T+hmax);
end

% Si tienen más filas (a veces guardas toda la historia), toma solo las primeras T+hmax o las últimas:
all_forecasts_cond = all_forecasts_cond(1:T+hmax, :, :);
all_forecasts_unc  = all_forecasts_unc(1:T+hmax, :, :);


valid_draws_cond = ~all(isnan(all_forecasts_cond(T+1,:,:)), [1 2]);
valid_draws_unc  = ~all(isnan(all_forecasts_unc(T+1,:,:)),  [1 2]);

all_forecasts_cond = all_forecasts_cond(:,:,valid_draws_cond);
all_forecasts_unc  = all_forecasts_unc(:,:,valid_draws_unc);

num_valid_cond = size(all_forecasts_cond,3);
num_valid_unc  = size(all_forecasts_unc,3);

if num_valid_cond == 0 || num_valid_unc == 0
    error('No valid forecast draws found. (Cond valid: %d, Unc valid: %d).', num_valid_cond, num_valid_unc);
end

cond_median = median(all_forecasts_cond, 3, 'omitnan');
cond_p16    = prctile(all_forecasts_cond, 16, 3);
cond_p84    = prctile(all_forecasts_cond, 84, 3);
unc_median  = median(all_forecasts_unc,  3, 'omitnan');

%% --- Ejes JOINED (incluyen T) ---
% dates_joined_model = [dates_historical(T); dates_forecast(:)];  % (hmax+1)x1
% 
% if use_xls
%     dates_joined_xls = [dates_historical(T); dates_fcast_xls(:)]; % (horizon_fcast+1)x1
% 
%     % Intersección FECHAS => “menor común”
%     [dates_common, ia_model, ia_xls] = intersect(dates_joined_model, dates_joined_xls, 'stable');
% else
%     dates_common = dates_joined_model;
%     ia_model = (1:length(dates_joined_model))';
%     ia_xls = [];
% end
dates_joined_model = [dates_historical(T); dates_forecast(:)];  % (hmax+1)

if use_xls
    dates_joined_xls = [dates_historical(T); dates_fcast_xls(:)]; % (horizon_fcast+1)

    % Menor común por HORIZONTE (evita problemas de stamps de fecha)
    H_common = 1 + min(hmax, length(dates_fcast_xls));   % incluye T

    dates_common = dates_joined_model(1:H_common);

    ia_model = (1:H_common)';
    ia_xls   = (1:H_common)';  % alineación “por paso” (T, T+1, ...)
else
    dates_common = dates_joined_model;
    ia_model = (1:length(dates_joined_model))';
    ia_xls = [];
end


%% --- Plotting ---
fig_handles = [];
max_per_fig = 9;
ncols = 3; nrows = 3;

for i = 1:nvar
    if mod(i-1, max_per_fig) == 0
        h_fig = figure('Name', sprintf('Conditional vs Unconditional Forecasts (Page %d)', ceil(i / max_per_fig)));
        fig_handles = [fig_handles; h_fig];
    end

    subplot_idx = mod(i-1, max_per_fig) + 1;
    subplot(nrows, ncols, subplot_idx);
    hold on;

    % ===== 1) Construir series joined del MODELO y recortar al común =====
    cond_med_joined = cond_median(T:end, i);
    cond_p16_joined = cond_p16(T:end, i);
    cond_p84_joined = cond_p84(T:end, i);
    unc_med_joined  = unc_median(T:end, i);

    cond_med_c = cond_med_joined(ia_model);
    cond_p16_c = cond_p16_joined(ia_model);
    cond_p84_c = cond_p84_joined(ia_model);
    unc_med_c  = unc_med_joined(ia_model);

    % % ===== 2) Banda 68% condicional =====
    fill([dates_common; flipud(dates_common)], ...
         [cond_p16_c; flipud(cond_p84_c)], ...
         [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % ===== 3) Historia =====
    plot(time_hist_win, yt_win(:, i), 'k-', 'LineWidth', 1.5);

    % ===== 4) Mediana condicional =====
    plot(dates_common, cond_med_c, 'b-', 'LineWidth', 1.5);

    % ===== 5) Mediana incondicional (del run condicional) =====
    plot(dates_common, unc_med_c, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

    % ===== 6) Excel Forecast_draws (banda + mediana) recortado al común =====
    if use_xls
        xls_med_joined = [yt(T,i); F_med_xls(:,i)];
        xls_p16_joined = [yt(T,i); F_p16_xls(:,i)];
        xls_p84_joined = [yt(T,i); F_p84_xls(:,i)];

        xls_med_c = xls_med_joined(ia_xls);
        xls_p16_c = xls_p16_joined(ia_xls);
        xls_p84_c = xls_p84_joined(ia_xls);

        % fill([dates_common; flipud(dates_common)], ...
        %      [xls_p16_c; flipud(xls_p84_c)], ...
        %      [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.35);

        plot(dates_common, xls_med_c, '--', 'LineWidth', 1.5);
    end

    % ===== 7) Condition path (alineado al común) =====
    if ~isempty(observableCondicionado_idx) && i == observableCondicionado_idx && ~all(isnan(yCondition(:, i)))
        dates_condpath_joined = [dates_historical(T); dates_forecast(:)];
        condpath_joined = [yt(T,i); yCondition(:, i)];

        [~, ia_cp, ia_dc] = intersect(dates_common, dates_condpath_joined, 'stable');
        if ~isempty(ia_cp)
            plot(dates_common(ia_cp), condpath_joined(ia_dc), 'g--', 'LineWidth', 1.5);
        end
    end

    % ===== 8) Realized (alineado al común) =====
    if show_realized && ~isempty(y_realized_future) && ~isempty(dates_realized_future)
        dates_real_joined = [dates_historical(T); dates_realized_future(:)];
        y_real_joined     = [yt(T,i); y_realized_future(:, i)];

        [~, ia_r, ia_dr] = intersect(dates_common, dates_real_joined, 'stable');
        if ~isempty(ia_r)
            plot(dates_common(ia_r), y_real_joined(ia_dr), 'r-.', 'LineWidth', 1.5);
        end
    end

    hold off;

    title(VarNames{i});
    xlabel('Fecha'); ylabel('Valor');
    grid on;

    % Eje X: historia acotada hasta fin del común
    % xlim([time_hist_win(1) dates_common(end)]);
    % datetick('x', 'yyyy-qq', 'keepticks');
    xlim([time_hist_win(1) dates_common(end)]);

    % Ticks trimestrales desde el inicio visible hasta el final común
    tick_start = dateshift(time_hist_win(1), 'start', 'quarter');
    tick_end   = dateshift(dates_common(end), 'start', 'quarter');
    xt = (tick_start : calquarters(1) : tick_end)';
    
    xticks(xt);
    xticklabels(compose('%d-Q%d', year(xt), quarter(xt)));

end

%% --- Leyenda maestra ---
try
    master_labels = {'Cond. 68% CI', 'History', 'Cond. Median', 'Uncond. Median'};
    if use_xls
        master_labels = [master_labels, {'ForecastDraws 68% CI', 'ForecastDraws Median'}];
    end
    master_labels = [master_labels, {'Condition Path'}];

    if show_realized
        master_labels = [master_labels, {'Realized'}];
    end

    hL = axes('Visible','off', 'Position', [0.15 0.1 0.7 0.05], 'Parent', fig_handles(end));
    hold(hL, 'on');

    dummy_plots = gobjects(1, length(master_labels));
    k = 1;

    dummy_plots(k) = fill(hL, [NaN NaN NaN NaN], [NaN NaN NaN NaN], [0.8 0.8 1], 'EdgeColor','none'); k=k+1;
    dummy_plots(k) = plot(hL, NaN, NaN, 'k-', 'LineWidth', 1.5); k=k+1;
    dummy_plots(k) = plot(hL, NaN, NaN, 'b-', 'LineWidth', 1.5); k=k+1;
    dummy_plots(k) = plot(hL, NaN, NaN, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); k=k+1;

    if use_xls
        dummy_plots(k) = fill(hL, [NaN NaN NaN NaN], [NaN NaN NaN NaN], [0.9 0.9 0.9], 'EdgeColor','none'); k=k+1;
        dummy_plots(k) = plot(hL, NaN, NaN, '--', 'LineWidth', 1.5); k=k+1;
    end

    dummy_plots(k) = plot(hL, NaN, NaN, 'g--', 'LineWidth', 1.5); k=k+1;

    if show_realized
        dummy_plots(k) = plot(hL, NaN, NaN, 'r-.', 'LineWidth', 1.5);
    end

    legend(dummy_plots, master_labels, 'Orientation','horizontal', 'Location','north', 'Box','off');
    hold(hL, 'off');
catch ME
    warning('No se pudo crear la leyenda única de pronósticos. Error: %s', ME.message);
end

% sgtitle('Conditional vs Unconditional Forecast (Median & 68% CI) + ForecastDraws');

end
