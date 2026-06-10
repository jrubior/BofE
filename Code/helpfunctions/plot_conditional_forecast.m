function [fig_handles] = plot_conditional_forecast(yf_save, yf_unc_save, yt, T, hmax, nvar, yCondition, VarNames, y_realized_future, dates_realized_future, dates_historical, dates_forecast, observableCondicionado_idx)
% PLOT_CONDITIONAL_FORECAST Grafica pronósticos condicionales vs. incondicionales usando fechas.
%
% VERSIÓN FINAL (13-Nov-2025):
%   - DEVUELVE [fig_handles] para un guardado fiable.
%   - MANTIENE el "salto" visual (no conecta historia con pronóstico).
%   - Incluye fechas (datetime) en el eje X.
%   - Incluye datos "Realized" como línea roja.
%   - Corrige el orden xlim / datetick.
%
% INPUTS:
%   yf_save:            Cell array con draws del pronóstico condicional.
%   yf_unc_save:        Cell array con draws del pronóstico incondicional.
%   yt:                 Datos históricos (T x nvar) usados para estimar.
%   T:                  Número de observaciones en yt (ya descuenta nlag).
%   hmax:               Horizonte del pronóstico.
%   nvar:               Número de variables.
%   yCondition:         Matriz (hmax x nvar) de condiciones sobre observables.
%   VarNames:           Cell array (1 x nvar) con nombres de variables.
%   y_realized_future:  Datos futuros observados (para comparación).
%   dates_realized_future: Vector datetime para y_realized_future.
%   dates_historical:   Vector datetime (T x 1) para el eje X de 'yt'.
%   dates_forecast:     Vector datetime (hmax x 1) para el eje X del pronóstico.
%
% OUTPUTS:
%   fig_handles:        Array de handles de las figuras creadas.
%

% --- Validar que los vectores de fechas coincidan con T y hmax ---
if length(dates_historical) ~= T
    error('La longitud de dates_historical (%d) no coincide con T (%d).', length(dates_historical), T);
end
if length(dates_forecast) ~= hmax
    error('La longitud de dates_forecast (%d) no coincide con hmax (%d).', length(dates_forecast), hmax);
end

num_draws = length(yf_save);

% --- Mostrar solo los últimos 15 períodos históricos ---
hist_len = 15;
hist_idx = max(1, T - hist_len + 1) : T;
time_hist_win = dates_historical(hist_idx); % Eje X de fechas acotado
yt_win = yt(hist_idx, :);                   % Historia acotada

% --- Process CONDITIONAL Forecasts ---
all_forecasts_cond = nan(T + hmax, nvar, num_draws);
empty_cond_draws = 0;
for d = 1:num_draws
    if ~isempty(yf_save{d}) && size(yf_save{d}, 1) == T + hmax && size(yf_save{d}, 2) == nvar
        all_forecasts_cond(:,:,d) = yf_save{d};
    else
        empty_cond_draws = empty_cond_draws + 1;
    end
end
valid_draws_cond = ~all(isnan(all_forecasts_cond(T+1,:,:)), [1 2]);
all_forecasts_cond = all_forecasts_cond(:,:,valid_draws_cond);
num_valid_cond = size(all_forecasts_cond, 3);

% --- Process UNCONDITIONAL Forecasts ---
all_forecasts_unc = nan(T + hmax, nvar, num_draws);
empty_unc_draws = 0;
for d = 1:num_draws
    if ~isempty(yf_unc_save{d}) && size(yf_unc_save{d}, 1) == T + hmax && size(yf_unc_save{d}, 2) == nvar
        all_forecasts_unc(:,:,d) = yf_unc_save{d};
    else
        empty_unc_draws = empty_unc_draws + 1;
    end
end
valid_draws_unc = ~all(isnan(all_forecasts_unc(T+1,:,:)), [1 2]);
all_forecasts_unc = all_forecasts_unc(:,:,valid_draws_unc);
num_valid_unc = size(all_forecasts_unc, 3);

if num_valid_cond == 0 || num_valid_unc == 0
    error('No valid forecast draws found. (Cond: %d, Unc: %d empty).', empty_cond_draws, empty_unc_draws);
end
fprintf('Loaded and processed %d valid conditional and %d valid unconditional forecast draws.\n', num_valid_cond, num_valid_unc);

% --- Calculate Percentiles ---
% Conditional
cond_median = median(all_forecasts_cond, 3, 'omitnan');
cond_p16 = prctile(all_forecasts_cond, 16, 3);
cond_p84 = prctile(all_forecasts_cond, 84, 3);
% Unconditional
unc_median = median(all_forecasts_unc, 3, 'omitnan');

% --- Plotting ---

% --- CAMBIO CLAVE: Inicializar el array de handles ---
fig_handles = []; 

max_per_fig = 9;        % máximo de subplots por figura
ncols = 3;              % 3 columnas
nrows = 3;              % 3 filas
var_names = VarNames;

for i = 1:nvar
     if mod(i-1, max_per_fig) == 0
        % --- CAMBIO CLAVE: Capturar el handle (h_fig) ---
        h_fig = figure('Name', sprintf('Conditional vs Unconditional Forecasts (Page %d)', ceil(i / max_per_fig)));
        fig_handles = [fig_handles; h_fig]; % Añadir el handle al array
    end
    
    % índice dentro de la figura actual (1–9)
    subplot_idx = mod(i-1, max_per_fig) + 1;
    subplot(nrows, ncols, subplot_idx);
    hold on;
    
    % 1. Plot CONDITIONAL credible interval bands (usando dates_forecast)
    % (Inicia en T+1)
    fill([dates_forecast(:); flipud(dates_forecast(:))], ...
         [cond_p16(T+1:end, i); flipud(cond_p84(T+1:end, i))], ...
         [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
         
    % 2. Plot historical data (usando time_hist_win)
    plot(time_hist_win, yt_win(:, i), 'k-', 'LineWidth', 1.5);
    
    % 3. Plot CONDITIONAL median forecast (usando dates_forecast)
    % (Inicia en T+1)
    plot(dates_forecast, cond_median(T+1:end, i), 'b-', 'LineWidth', 1.5); % Blue
    
    % 4. Plot UNCONDITIONAL median forecast (usando dates_forecast)
    % (Inicia en T+1)
    plot(dates_forecast, unc_median(T+1:end, i), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); % Gray dotted
    
    % 5. Plot the conditioning path (if applicable)
    if ~isempty(observableCondicionado_idx) && i == observableCondicionado_idx && ~all(isnan(yCondition(:, i)))
        plot(dates_forecast, yCondition(:, i), 'g--', 'LineWidth', 1.5); % Red dashed
    end
    
    % 6. (Opcional) Plotear datos futuros realizados (línea roja)
    if ~isempty(y_realized_future) && ~isempty(dates_realized_future)
        plot_realized = y_realized_future(:, i);
        plot_dates = dates_realized_future;
        idx_realized = (plot_dates >= dates_forecast(1)) & (plot_dates <= dates_forecast(end));
        if any(idx_realized)
            plot(plot_dates(idx_realized), plot_realized(idx_realized), 'r-.', 'LineWidth', 1.5); % Línea roja
        end
    end

    hold off;
    
    title(var_names{i});
    
    xlabel('Fecha');
    ylabel('Value');
    grid on;
    
    % 1. Asegurar que el eje X cubra el historial y el pronóstico PRIMERO
    if ~isempty(time_hist_win)
        xlim([time_hist_win(1) dates_forecast(end)]);
    else
        xlim([dates_forecast(1) dates_forecast(end)]);
    end
    
    % 2. AHORA, formatear el eje X
    datetick('x', 'yyyy-qq', 'keepticks'); 
end

% --- Creación de leyenda única "Maestra" para toda la figura (CORREGIDO) ---
try
    % 1. Definir la "leyenda maestra" con TODOS los elementos posibles
    master_labels = {'Cond. 68% CI', 'History', 'Cond. Median', 'Uncond. Median', 'Condition Path', 'Realized'};
    
    % 2. Crear un eje invisible en la parte superior
    hL = axes('Visible','off', 'Position', [0.15 0.1 0.7 0.05], 'Parent', fig_handles(end));
    hold(hL, 'on');
    
    % 3. Dibujar líneas "dummy" Y CAPTURAR SUS HANDLES en un array
    %    (Ya no usamos 'HandleVisibility', 'off')
    dummy_plots = gobjects(1, 6); % Pre-alocar array de handles
    dummy_plots(1) = fill(hL, [NaN NaN NaN NaN], [NaN NaN NaN NaN], [0.8 0.8 1], 'EdgeColor', 'none');
    dummy_plots(2) = plot(hL, NaN, NaN, 'k-', 'LineWidth', 1.5);
    dummy_plots(3) = plot(hL, NaN, NaN, 'b-', 'LineWidth', 1.5);
    dummy_plots(4) = plot(hL, NaN, NaN, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    dummy_plots(5) = plot(hL, NaN, NaN, 'g--', 'LineWidth', 1.5);
    dummy_plots(6) = plot(hL, NaN, NaN, 'r-.', 'LineWidth', 1.5);

    % 4. Crear la leyenda pasando AMBOS, los handles y las etiquetas
    %    Ahora los números coinciden (6 handles, 6 etiquetas)
    legend(dummy_plots, master_labels, ...
           'Orientation', 'horizontal', 'Location', 'north', 'Box', 'off');
    hold(hL, 'off');
catch ME
    warning('No se pudo crear la leyenda única de pronósticos. Error: %s', ME.message);
end
% ===================================

sgtitle('Conditional vs Unconditional Forecast (Median & 68% CI)');

end % Fin de la función