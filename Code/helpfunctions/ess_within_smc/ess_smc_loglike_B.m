function loglik_prop = ess_smc_loglike_B(z_prop,function_restrictions_Svec,cutoff,Sigmadraw,Qdraw,info)
% likelihood eval y given z_prop

% --- recover B
B_prop = reshape(z_prop,info.m,info.nvar);

% --- Check restrictions
Svec_prop = function_restrictions_Svec(B_prop, Sigmadraw, Qdraw, info);
S_prop    = all(cutoff < Svec_prop);

% --- Likelihood 
loglik_prop = -Inf;
if S_prop==1
    loglik_prop = log(1);%logpdfBgivenSigma ;
end
