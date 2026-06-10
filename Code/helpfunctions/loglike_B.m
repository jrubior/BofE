function loglik_prop = loglike_B(z_prop,function_restrictions,fo_inv,fo_str2irfs,Sigmadraw,Qdraw,info,mmuTilde,OomegaTilde)
% likelihood eval y given z_prop
% z_prop = [ddelta, ggamma]


% ---

B_prop = reshape(z_prop,info.m,info.nvar);




S_prop = function_restrictions(B_prop, Sigmadraw, Qdraw,info);




% p(B|Sigma)
%varBgivenSigma = kron(Sigmadraw,OomegaTilde);
%varBgivenSigma = (varBgivenSigma+varBgivenSigma')/2;
 
%logpdfBgivenSigma = logmvnpdf_mc(z_prop,vec(mmuTilde),varBgivenSigma);%log(mvnpdf(vec(Bdraw),vec(mmuTilde),varBgivenSigma));


loglik_prop = -Inf;
if S_prop==1
loglik_prop = log(1);%logpdfBgivenSigma ;
end
