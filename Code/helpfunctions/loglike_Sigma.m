function loglik_prop = loglike_Sigma(z_prop,function_restrictions,fo_inv,fo_str2irfs,Bdraw,Qdraw,info,nnuTilde,mmuTilde,OomegaTilde)
% likelihood eval y given z_prop
% z_prop = [ddelta, ggamma]


% ---

R_prop = reshape(z_prop,info.nvar,nnuTilde);


Sigma_prop=inv(R_prop*R_prop');


S_prop = function_restrictions(Bdraw, Sigma_prop, Qdraw,fo_inv,fo_str2irfs,info);




% p(B|Sigma)
varBgivenSigma = kron(Sigma_prop,OomegaTilde);
% varBgivenSigma = (varBgivenSigma+varBgivenSigma')/2;
 

logpdfBgivenSigma = logmvnpdf_mc(vec(Bdraw),vec(mmuTilde),varBgivenSigma);%log(mvnpdf(vec(Bdraw),vec(mmuTilde),varBgivenSigma));

loglik_prop = -Inf;
if S_prop==1
    
loglik_prop = logpdfBgivenSigma ;
end
