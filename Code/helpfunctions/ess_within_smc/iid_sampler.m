function [Bdraw, Sigmadraw, Rdraw_vec, Qdraw, Xdraw_vec, cholSigmadraw] = iid_sampler(mmuTilde, cholOomegaTilde, inv_cholPpsiTilde_prime, nnuTilde, info)
% --- generate B,Sigma,Q from the posterior distribution
% Sigmadraw     = iwishrnd(PpsiTilde,nnuTilde);

% Sigma and R
Rdraw = inv_cholPpsiTilde_prime * randn(info.nvar, nnuTilde);
Sigmadraw = eye(info.nvar) / (Rdraw * Rdraw');
Rdraw_vec = Rdraw(:);

% Bdraw
% cholSigmadraw = chol(Sigmadraw)';
% Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(info.m*info.nvar,1) + reshape(mmuTilde,info.nvar*info.m,1);
% Bdraw         = reshape(Bdraw,info.nvar*info.nlag+info.nex,info.nvar);
Bdraw = cholOomegaTilde*randn(info.pred,info.nvar)*chol(Sigmadraw) + mmuTilde;

% Q and X
[Qdraw,~,Xdraw]  = DrawQ(info.nvar);
Xdraw_vec = Xdraw(:);