function loglik_prop = ess_smc_loglike_Q(z_prop, function_restrictions_Svec, cutoff, Bdraw, Sigmadraw, info)
% likelihood eval y given z_prop

% --- recover Q
X_prop = reshape(z_prop,info.nvar,info.nvar);
[Q_prop,~] = qr_unique(X_prop);

% --- Check restrictions
Svec_prop = function_restrictions_Svec(Bdraw, Sigmadraw, Q_prop, info);
S_prop    = all(cutoff < Svec_prop);

% --- Likelihood 
loglik_prop = -Inf;
if S_prop==1 
    loglik_prop = log(1);
end