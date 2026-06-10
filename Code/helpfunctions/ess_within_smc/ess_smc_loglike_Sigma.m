function loglik_prop = ess_smc_loglike_Sigma(z_prop,function_restrictions_Svec,cutoff,Bdraw,Qdraw,info,nnuTilde,mmuTilde,chol_lower_OomegaTilde,chol_lower_OomegaTilde_inv)
% likelihood eval y given z_prop
% z_prop = [ddelta, ggamma]


% ---
R_prop = reshape(z_prop,info.nvar,nnuTilde);
Sigma_prop=eye(info.nvar)/(R_prop*R_prop');

% Check restriction
Svec_prop = function_restrictions_Svec(Bdraw, Sigma_prop, Qdraw, info);
S_prop    = all(cutoff < Svec_prop);


% --- Likelihood calculation when restrictions hold
loglik_prop = -Inf;
if S_prop==1

    % ---
    % p(B|Sigma)
    % ---
    % Fast likelihood calculation using kronecker product structure
    % note that
    % chol_upper(V) = kron(chol_upper(Sigma),chol_upper(Omega))
    % Want to find xRinv such that  xRinv * chol_upper(V) = X0 where xRinv
    % and X0 are row vectors
    % Then, chol_lower(V) * xRinv' = X0'
    % Then, chol_lower(Omega) * unvec(xRinv') * chol_lower(Sigma)' = unvec(X0')
    % And therefore, unvec(xRinv') = chol_lower(Omega) \ unvec(X0') / chol_lower(Sigma)'
    % dim(Omega) = n_Omega (info.m)
    % dim(Sigma) = n_Sigma (info.nvar)
    % ---
    chol_lower_Sigma_prop = chol(Sigma_prop,'lower');
    % 
    % % actual calculation
    % unvec_X0_prime = reshape( ( vec(Bdraw) - vec(mmuTilde) )', info.m, info.nvar);
    % % unvec_xRinv_prime = chol_lower_OomegaTilde \ unvec_X0_prime / chol_lower_Sigma_prop;
    % unvec_xRinv_prime = chol_lower_OomegaTilde_inv * unvec_X0_prime / chol_lower_Sigma_prop;
    % xRinv_prime = unvec_xRinv_prime(:);
    % 
    % logSqrtDetSigma = sum(log(diag(chol_lower_Sigma_prop))) * info.nvar + sum(log(diag(chol_lower_OomegaTilde))) * info.m;
    % 
    % quadform = sum(xRinv_prime.^2, 1); %xRinv_prime is column vector
    % 
    % logpdfBgivenSigma = (-0.5*quadform - logSqrtDetSigma - size(xRinv_prime,1)*log(2*pi)/2);

    % % actual calculation: jonas
unvec_X0_prime = reshape( ( vec(Bdraw) - vec(mmuTilde) ), info.m, info.nvar);
% unvec_xRinv_prime = chol_lower_OomegaTilde \ unvec_X0_prime / chol_lower_Sigma_prop;
unvec_xRinv_prime = chol_lower_OomegaTilde_inv' * unvec_X0_prime / chol_lower_Sigma_prop;
xRinv_prime = unvec_xRinv_prime(:);

logSqrtDetSigma = sum(log(diag(chol_lower_Sigma_prop))) * info.m + sum(log(diag(chol_lower_OomegaTilde))) * info.nvar;

quadform = sum(xRinv_prime.^2, 1); %xRinv_prime is column vector

logpdfBgivenSigma = (-0.5*quadform - 0.5*logSqrtDetSigma - size(xRinv_prime,1)*log(2*pi)/2);



    loglik_prop = logpdfBgivenSigma ;
end
