Sigma = [2 0.5; 0.5 1.5];      % 2x2 positive definite
Omega = [1.2 0.3; 0.3 0.8];    % 2x2 positive definite

info.nvar = 2;   % dimension of Sigma
info.m = 2;      % dimension of Omega

Bdraw = [0.5 1.0; -0.5 0.2];    % 2x2 matrix
mmuTilde = [0.4 0.9; -0.4 0.1]; % prior mean (2x2 matrix)

% Vectorize
vec_Bdraw = Bdraw(:);
vec_mmuTilde = mmuTilde(:);

% Cholesky factors
chol_lower_Sigma = chol(Sigma, 'lower');
chol_lower_Omega = chol(Omega, 'lower');

% Calculate unvec(X0') for your code
unvec_X0_prime = reshape((vec_Bdraw - vec_mmuTilde)', info.m, info.nvar);

% Solve for unvec(xRinv')
unvec_xRinv_prime = chol_lower_Omega \ unvec_X0_prime / chol_lower_Sigma;

% Reshape to vector
xRinv_prime = unvec_xRinv_prime(:);

% Calculate logSqrtDetSigma using your formula
logSqrtDetSigma = sum(log(diag(chol_lower_Sigma))) * info.m + sum(log(diag(chol_lower_Omega))) * info.nvar;

% Calculate quadform
quadform = sum(xRinv_prime.^2);

% Calculate log-likelihood
logpdfBgivenSigma = (-0.5*quadform - logSqrtDetSigma - size(xRinv_prime,1)*log(2*pi)/2);


V = kron(Sigma, Omega);

% Standard full multivariate normal log-likelihood
diff = vec_Bdraw - vec_mmuTilde;
quadform_naive = diff' * (V \ diff);

logdetV = log(det(V));
d = length(vec_Bdraw);

logpdf_naive = -0.5 * (quadform_naive + logdetV + d*log(2*pi));

% maltab
logpdf_matlab = logmvnpdf_mc(vec_Bdraw, vec_mmuTilde,V);

% 
disp(['Your code logpdf = ', num2str(logpdfBgivenSigma)]);
disp(['Naive full V logpdf = ', num2str(logpdf_naive)]);
disp(['Naive full V logpdf = ', num2str(logpdf_matlab)]);

%%
Sigma = 2; 
Omega = 3;

info.nvar = 1;
info.m = 1;

Bdraw = 1.5;
mmuTilde = 1.0;

vec_Bdraw = Bdraw(:);
vec_mmuTilde = mmuTilde(:);

chol_lower_Sigma = chol(Sigma, 'lower');
chol_lower_Omega = chol(Omega, 'lower');

unvec_X0_prime = (vec_Bdraw - vec_mmuTilde)'; % row vector
unvec_xRinv_prime = chol_lower_Omega \ unvec_X0_prime / chol_lower_Sigma;
xRinv_prime = unvec_xRinv_prime(:);

logSqrtDetSigma = sum(log(diag(chol_lower_Sigma))) * info.m + sum(log(diag(chol_lower_Omega))) * info.nvar;

quadform = sum(xRinv_prime.^2);

logpdfBgivenSigma = (-0.5*quadform - logSqrtDetSigma - size(xRinv_prime,1)*log(2*pi)/2);

disp(['Logpdf = ', num2str(logpdfBgivenSigma)]);

