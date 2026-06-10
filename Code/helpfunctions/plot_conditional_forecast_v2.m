function [fig_handles] = plot_conditional_forecast_v2(yf_save, yf_unc_save, yt, T, hmax, nvar, yCondition, VarNames, y_realized_future, dates_realized_future, dates_historical, dates_forecast, observableCondicionado_idx, show_realized, num_hist_periods)
% PLOT_CONDITIONAL_FORECAST_V2 Grafica pronósticos condicionales vs. incondicionales con mayor flexibilidad.
%
% CAMBIOS SOLICITADOS:
% 1. Control flexible para mostrar/ocultar 'Realized' (show_realized).
% 2. Control sobre el número de períodos históricos a mostrar (num_hist_periods).
% 3. Unir visualmente la historia con el pronóstico (añadiendo T a la línea de pronóstico).
%
% INPUTS ADICIONALES:
%   show_realized:      Booleano (true/false) para plotear y_realized_future.
%   num_hist_periods:   Número de períodos históricos a mostrar.
%
% (Resto de INPUTS son los mismos que en la versión original)
%
% OUTPUTS:
%   fig_handles:        Array de handles de las figuras creadas.

% --- Validar que los vectores de fechas coincidan con T y hmax ---
if length(dates_historical) ~= T
    error('La longitud de dates_historical (%d) no coincide con T (%d).', length(dates_historical), T);
end
if length(dates_forecast) ~= hmax
    error('La longitud de dates_forecast (%d) no coincide con hmax (%d).', length(dates_forecast), hmax);
end

% 2. Control sobre la historia mostrada (Acortar los datos reales)
%    Usamos el nuevo argumento num_hist_periods
hist_len = num_hist_periods; 
hist_idx = max(1, T - hist_len + 1) : T;
time_hist_win = dates_historical(hist_idx);  % Eje X de fechas acotado
yt_win = yt(hist_idx, :);                    % Historia acotada

% --- Process Forecasts (Mismo proceso de carga y validación) ---
num_draws = length(yf_save);
all_forecasts_cond = nan(T + hmax, nvar, num_draws);
all_forecasts_unc = nan(T + hmax, nvar, num_draws);

% Carga Condicional
for d = 1:num_draws
    if ~isempty(yf_save{d}) && size(yf_save{d}, 1) == T + hmax && size(yf_save{d}, 2) == nvar
        all_forecasts_cond(:,:,d) = yf_save{d};
    end
end
valid_draws_cond = ~all(isnan(all_forecasts_cond(T+1,:,:)), [1 2]);
all_forecasts_cond = all_forecasts_cond(:,:,valid_draws_cond);
num_valid_cond = size(all_forecasts_cond, 3);

% Carga Incondicional
for d = 1:num_draws
    if ~isempty(yf_unc_save{d}) && size(yf_unc_save{d}, 1) == T + hmax && size(yf_unc_save{d}, 2) == nvar
        all_forecasts_unc(:,:,d) = yf_unc_save{d};
    end
end
valid_draws_unc = ~all(isnan(all_forecasts_unc(T+1,:,:)), [1 2]);
all_forecasts_unc = all_forecasts_unc(:,:,valid_draws_unc);
num_valid_unc = size(all_forecasts_unc, 3);

if num_valid_cond == 0 || num_valid_unc == 0
    error('No valid forecast draws found. (Cond: %d, Unc: %d).', num_draws - num_valid_cond, num_draws - num_valid_unc);
end
fprintf('Loaded and processed %d valid conditional and %d valid unconditional forecast draws.\n', num_valid_cond, num_valid_unc);

% --- Calculate Percentiles (Mismo cálculo) ---
cond_median = median(all_forecasts_cond, 3, 'omitnan');
cond_p16 = prctile(all_forecasts_cond, 16, 3);
cond_p84 = prctile(all_forecasts_cond, 84, 3);
unc_median = median(all_forecasts_unc, 3, 'omitnan');

% --- Preparación para Unir Historia y Pronóstico (Punto 3) ---
% Fechas: Incluimos el último punto de la historia (T) y luego el forecast (T+1 a T+hmax)
% Línea 71 Corregida:
dates_forecast_joined = [dates_historical(T); dates_forecast']; % (hmax + 1) x 1
% Datos: Incluimos la observación histórica T al inicio del forecast
% Los arrays de medianas ya contienen la observación T en la posición T.
% Necesitamos cond_median(T:end, i) -> (hmax + 1) x 1
% Necesitamos cond_pXX(T:end, i)   -> (hmax + 1) x 1

% --- Plotting ---
fig_handles = []; 
max_per_fig = 9;        
ncols = 3;              
nrows = 3;              
var_names = VarNames;

for i = 1:nvar
    if mod(i-1, max_per_fig) == 0
        h_fig = figure('Name', sprintf('Conditional vs Unconditional Forecasts (Page %d)', ceil(i / max_per_fig)));
        fig_handles = [fig_handles; h_fig];
    end
    
    subplot_idx = mod(i-1, max_per_fig) + 1;
    subplot(nrows, ncols, subplot_idx);
    hold on;
    
    % 1. Plot CONDITIONAL credible interval bands (usando dates_forecast)
    %    Ahora va de T hasta T+hmax, y usamos dates_forecast_joined
    fill([dates_forecast_joined(:); flipud(dates_forecast_joined(:))], ...
         [cond_p16(T:end, i); flipud(cond_p84(T:end, i))], ...
         [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
         
    % 2. Plot historical data (usando time_hist_win)
    plot(time_hist_win, yt_win(:, i), 'k-', 'LineWidth', 1.5);
    
    % 3. Plot CONDITIONAL median forecast (UNIDO)
    %    Va de T hasta T+hmax. Usamos cond_median(T:end, i) y dates_forecast_joined.
    plot(dates_forecast_joined, cond_median(T:end, i), 'b-', 'LineWidth', 1.5); % Blue
    
    % 4. Plot UNCONDITIONAL median forecast (UNIDO)
    %    Va de T hasta T+hmax. Usamos unc_median(T:end, i) y dates_forecast_joined.
    plot(dates_forecast_joined, unc_median(T:end, i), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); % Gray dotted
    
    % 5. Plot the conditioning path (if applicable)
    %    NOTA: La condición es solo para T+1 a T+hmax, NO T. Usamos dates_forecast.
    if ~isempty(observableCondicionado_idx) && i == observableCondicionado_idx && ~all(isnan(yCondition(:, i)))
        plot(dates_forecast, yCondition(:, i), 'g--', 'LineWidth', 1.5); % Green dashed
    end
    
    % --- 6. (Opcional) Plotear datos futuros realizados (Punto 1: Control)
if show_realized && ~isempty(y_realized_future) && ~isempty(dates_realized_future)
    
    % 1. OBTENER EL PUNTO DE CONEXIÓN HISTÓRICO (Última observación de estimación)
    last_historical_value = yt(T, i);
    last_estimation_date = dates_historical(T);
    
    % 2. UNIR DATOS: last_historical_value + y_realized_future (continuidad)
    % Se usan los vectores columna (:) para la concatenación vertical segura
    y_realized_joined = [last_historical_value; y_realized_future(:, i)];
    dates_realized_joined = [last_estimation_date; dates_realized_future]; 
    
    plot_realized = y_realized_joined; % Usamos el vector unido
    plot_dates = dates_realized_joined; % Usamos las fechas unidas
    
    % Solo plotea el segmento que se superpone con el horizonte de pronóstico
    % La condición de índice ahora debe comenzar en el punto de conexión (last_estimation_date)
    idx_realized = (plot_dates >= dates_forecast_joined(1)) & (plot_dates <= dates_forecast_joined(end));
    
    if any(idx_realized)
        plot(plot_dates(idx_realized), plot_realized(idx_realized), 'r-.', 'LineWidth', 1.5); % Línea roja
    end
end
    hold off;
    
    title(var_names{i});
    
    xlabel('Fecha');
    ylabel('Valor');
    grid on;
    
    % 1. Asegurar que el eje X cubra la historia acotada y el pronóstico
    if ~isempty(time_hist_win)
        % El límite termina en el último punto del pronóstico
        xlim([time_hist_win(1) dates_forecast_joined(end)]); 
    else
        xlim([dates_forecast_joined(1) dates_forecast_joined(end)]);
    end
    
    % 2. Formatear el eje X
    datetick('x', 'yyyy-qq', 'keepticks'); 
end

% --- Creación de leyenda única "Maestra" para toda la figura ---
try
    % 1. Definir la "leyenda maestra" con TODOS los elementos posibles
    master_labels = {'Cond. 68% CI', 'History', 'Cond. Median', 'Uncond. Median', 'Condition Path'};
    
    % Ajustar la leyenda si se muestra la ruta 'Realized'
    if show_realized
        master_labels = [master_labels, {'Realized'}];
    end
    
    % 2. Crear un eje invisible en la última figura
    hL = axes('Visible','off', 'Position', [0.15 0.1 0.7 0.05], 'Parent', fig_handles(end));
    hold(hL, 'on');
    
    % 3. Dibujar líneas "dummy" y CAPTURAR SUS HANDLES
    dummy_plots = gobjects(1, length(master_labels)); 
    
    dummy_plots(1) = fill(hL, [NaN NaN NaN NaN], [NaN NaN NaN NaN], [0.8 0.8 1], 'EdgeColor', 'none'); % Cond. 68% CI
    dummy_plots(2) = plot(hL, NaN, NaN, 'k-', 'LineWidth', 1.5); % History
    dummy_plots(3) = plot(hL, NaN, NaN, 'b-', 'LineWidth', 1.5); % Cond. Median
    dummy_plots(4) = plot(hL, NaN, NaN, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); % Uncond. Median
    dummy_plots(5) = plot(hL, NaN, NaN, 'g--', 'LineWidth', 1.5); % Condition Path

    if show_realized
        dummy_plots(6) = plot(hL, NaN, NaN, 'r-.', 'LineWidth', 1.5); % Realized
    end

    % 4. Crear la leyenda
    legend(dummy_plots, master_labels, ...
           'Orientation', 'horizontal', 'Location', 'north', 'Box', 'off');
    hold(hL, 'off');
catch ME
    warning('No se pudo crear la leyenda única de pronósticos. Error: %s', ME.message);
end

sgtitle('Conditional vs Unconditional Forecast (Median & 68% CI)');

end % Fin de la función