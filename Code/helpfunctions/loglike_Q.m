function loglik_prop = loglike_Q(z_prop,function_restrictions,gs_qr,fo_inv,fo_str2irfs,Bdraw,Sigmadraw,info)
% likelihood eval y given z_prop
% z_prop = [ddelta, ggamma]


% ---

X_prop = reshape(z_prop,info.nvar,info.nvar);


[Q_prop,~]=gs_qr(X_prop);

S_prop = function_restrictions(Bdraw, Sigmadraw, Q_prop,info);


loglik_prop = -Inf;


if S_prop==1 
    loglik_prop = log(1);
end




