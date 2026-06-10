function [Bdraw,Sigmadraw,Rdraw,Qdraw,Xdraw,n_try_B, n_try_S, n_try_Q] = ess_sampler(Bdraw,Sigmadraw,Rdraw,Qdraw,Xdraw,function_restrictions_Svec, cutoff, info, info_ess)
% -------------------------------------------------------------------------
% ESS sampler within SMC
% -------------------------------------------------------------------------
% Goal: Given (B,Sigma,Q), it generates another (B,Sigma,Q) given restrictions
% using transition kernel defined by ESS
% -------------------------------------------------------------------------
% Inputs
% (Bdraw, Sigmadraw, Qdraw): current draw
% cutoff: current bounds for the sign restrictions (Restrictions: cutoff < Svec)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Draw Q
% -------------------------------------------------------------------------
slice.scale_z    = 1;
slice.mean       = zeros(info.nvar*info.nvar,1);
slice.chol_cov_z = 1;

slice.fcn_lik    = @(z_prop) ess_smc_loglike_Q(z_prop,function_restrictions_Svec,cutoff,Bdraw,Sigmadraw,info);

slice.nobs       = info.nvar*info.nvar;

% --- actual slice sampling
% z_old = randn(info.nvar*info.nvar,1);
lik_old = slice.fcn_lik(Xdraw);
[Xdraw, ~, n_try_Q] = slice_sampling_v02(slice, Xdraw, lik_old);
[Qdraw,~] =  qr_unique(reshape(Xdraw,info.nvar,info.nvar));


% -------------------------------------------------------------------------
% Draw Sigma
% -------------------------------------------------------------------------
slice.scale_z    = 1;
slice.mean       = vec(info_ess.Rmean);
slice.chol_cov_z_common = info_ess.chol_lower_Rvariance_common; %
slice.nobs1 = info.nvar;
slice.nobs2 = info_ess.nnuTilde;


% slice.fcn_lik    = @(z_prop) loglike_Sigma_mc(z_prop,function_restrictions_i,fo_inv,fo_str2irfs,Bdraw,Qdraw,info,nnuTilde,mmuTilde,chol_lower_OomegaTilde,chol_lower_OomegaTilde_inv);

slice.fcn_lik = @(z_prop) ess_smc_loglike_Sigma(z_prop,function_restrictions_Svec,cutoff,Bdraw,Qdraw,info,info_ess.nnuTilde,info_ess.mmuTilde,info_ess.chol_lower_OomegaTilde,info_ess.chol_lower_OomegaTilde_inv);

slice.nobs       = info.nvar*info_ess.nnuTilde;

% --- actual slice sampling
lik_old = slice.fcn_lik(Rdraw);

[Rdraw, ~, n_try_S] = slice_sampling_v02_Sigma(slice, Rdraw, lik_old);

Rdraw_mat = reshape(Rdraw,info.nvar,info_ess.nnuTilde);
Sigmadraw=inv(Rdraw_mat*Rdraw_mat');


% -------------------------------------------------------------------------
% Draw B
% -------------------------------------------------------------------------
cholLSigmadraw = chol(Sigmadraw,'lower');

slice.scale_z    = 1;
slice.mean       = vec(info_ess.mmuTilde);

slice.chol_cov_z1 = cholLSigmadraw;%chol(Rvariance,'lower');
slice.chol_cov_z2 = info_ess.cholOomegaTilde;%chol(Rvariance,'lower');
slice.nobs1 = info.nvar;
slice.nobs2 = info.m;

slice.fcn_lik    = @(z_prop) ess_smc_loglike_B(z_prop,function_restrictions_Svec,cutoff,Sigmadraw,Qdraw,info);
slice.nobs       = info.nvar*info.m;

% --- actual slice sampling
z_old_B = Bdraw(:);
lik_old = slice.fcn_lik(z_old_B);
[z_old_B, ~, n_try_B] = slice_sampling_v02_kron(slice, z_old_B, lik_old);
Bdraw = reshape(z_old_B,info.m,info.nvar);

%x            = [vec(Bdraw); vec(Sigmadraw); vec(Qdraw)];
% structpara   = f_h_inv(x,info);
% irfpara      = fo_str2irfs(structpara);


% if counter_Q==1000
%        counter_Q=0;
% end

