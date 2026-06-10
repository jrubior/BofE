function [opt_cutoff, opt_qt] = cutoff_search(mat_Svec, target_ar, eff_zero)
% mat_Svec: restriction vector (R by nsim) where R is number of restrictions
% target_ar: target acceptance rate
% eff_zero: effeictive zero (that is, small number that we regard it as zero)

% Goal: we search for the largest qt_cutoff that is larger than target_ar


% --- setup
ngrid = 100;

qt_grid = linspace(0.01,0.99,ngrid)';
ar_grid = zeros(ngrid,1);

% --- calculate ar per each quantile
for i=1:ngrid
    qt_cutoff = qt_grid(i);

    lv_cutoff = min( quantile(mat_Svec, qt_cutoff, 2), eff_zero); %actual cut-off (note we oriented restrictions so that 0 < Svec; if cutoff is larger than 0 then we let the cutoff to zero)
    ar_cutoff = mean(all(lv_cutoff < mat_Svec)); % proportion of survived draws (acceptance rate)
    ar_grid(i) = ar_cutoff;
end

% --- find the largest qt that leads to ar larger than target_ar
opt_ind = find(ar_grid > target_ar, 1, 'last');

opt_qt = qt_grid(opt_ind);
opt_cutoff = min( quantile(mat_Svec, opt_qt, 2), eff_zero); 
