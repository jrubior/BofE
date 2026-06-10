clearvars
close all

% test github

data = fullfile('data.csv');
opts = detectImportOptions(data);
opts = opts.setvartype("date", "datetime");
opts = opts.setvaropts("date", "InputFormat",'yyyy.Q');
tbl = readtimetable(data, opts);
Y = tbl{:,2:end}*100;


%% Model parameters 
nvar      = size(Y,2);       % number of endogenous variables
nlag      = 1;               % number of lags
nex       = 1;               % set equal to 1 if a constant is included; 0 otherwise
m         = nvar*nlag + nex; % number of exogenous variables
T         = size(Y,1) - nlag; % sample
pos       = []; % if there is variables that are stationary.

%% MN Dummy Prior Parameters

disp('MN Prior Dummies')

Vc=1e4;
lambda1=0.2;
lambda3=4;

% These are the lambdas as defined in dummies
llambda=zeros(3,1);
llambda(1,1)=lambda1;
llambda(2,1)=lambda3/2;
llambda(3,1)=1; % This needs to be == 1 if we want to have the same result as in MN as GLP

diagPphiBar=zeros(nvar,1);
for i=1:nvar
    Mdl = arima(nlag, 0, 0);       % AR(1) model: p=1, d=0, q=0
    EstMdl = estimate(Mdl, Y(:,i),'Display','off');
    diagPphiBar(i)=EstMdl.Variance;
end
pphi=diagPphiBar;

% model = struct('Vc',Vc,'m',m','nvar',nvar,'nlag',nlag,'llambda', llambda,'ppsi', ppsi);
%[nnuBar, PphiBar, PpsiBar, OomegaBar] = minnesota_prior_with_dummiesv2(model);

% Initialize matrices
Ystar = zeros(m + llambda(3,1)*nvar, nvar);
Xstar = zeros(m + llambda(3,1)*nvar, m);

% this is the posterior associate with this prior
% | \Sigma | ^{(v + nvar +1)/2}
% and Ystar and Xstar data

% Standard deviation vector
sBar = sqrt(pphi);

%% Vectorized construction of dummy observations
% First lag prior (first nvar rows)
% Create identity matrix scaled by (1/lambda1) * sBar
I_n = eye(nvar);
scaling_first_lag = (1/llambda(1,1)) * sBar;
Ystar(1:nvar, 1:nvar) = I_n .* scaling_first_lag;
Xstar(1:nvar, 1:nvar) = I_n .* scaling_first_lag;

% Other lags prior (rows nvar+1 to m-1, assuming constant is last)
for ell = 2:nlag
    row_idx = (ell-1)*nvar + (1:nvar);  % Vectorized row indices
    col_idx = (ell-1)*nvar + (1:nvar);  % Vectorized column indices

    % Decay factor for this lag
    decay_factor = ell^llambda(2,1);
    decay_factor = ell^(-llambda(2,1));
    scaling_other_lags = (1/llambda(1,1)) * sBar * decay_factor;

    Xstar(row_idx, col_idx) = I_n .* scaling_other_lags;
end

% Prior for the constant (assuming constant is the last variable)
Xstar(m, m) = 1/sqrt(Vc);

% Prior for covariance matrix Sigma (dummy observations)
for ii = 1:llambda(3,1)
    row_idx = m + (ii-1)*nvar + (1:nvar);  % Vectorized row indices
    Ystar(row_idx, 1:nvar) = I_n .* sBar;
end

% nnuBar
nnuBar = size(Ystar,1) - m +  2; % T + nnuBar of prior before dummies

% OomegaBar
OomegaBarInverse = Xstar' * Xstar;
OomegaBar = OomegaBarInverse \ eye(m);

% PpsiBar
PpsiBar = OomegaBar * (Xstar' * Ystar);

% PphiBar
PphiBar = Ystar' * Ystar - PpsiBar' * OomegaBarInverse * PpsiBar;
% Ensure PphiBar is symmetric
PphiBar = (PphiBar + PphiBar') / 2;


%% Number of Draws 

nd=1e4;

%% Model definition using conjugatebvarm

disp('conjugatebvarm')

PriorMdl = conjugatebvarm(...
    nvar, ...                                                    % Number of Series
    nlag, ...                                                    % Number of Lags
    SeriesNames     = tbl.Properties.VariableNames(2:end), ...   % Optional, name of the series
    Omega           = PphiBar, ...                               % Optional, Default is false
    DoF             = nnuBar, ...                                % PphiBar
    Mu              = PpsiBar, ...                               % PpsiBar
    V               = OomegaBar, ...                             % OomegaBar
    IncludeConstant = nex, ...                                   % Optional, Default is true
    IncludeTrend    = false);  

PriorMdl.Omega;
PriorMdl.DoF;

PostMdl = PriorMdl.estimate(Y, Display = 'off');

[Bdraws,Sigmadraws] = simulate(PriorMdl, Y, NumDraws = nd);

PostMdl.Omega;
PostMdl.DoF;

if T + PriorMdl.DoF == PostMdl.DoF
    fprintf('Degrees of freedom match: %d + %d = %d\n', T, PriorMdl.DoF, PostMdl.DoF);
else
    error('Degrees of freedom mismatch: %d + %d ≠ %d', T, PriorMdl.DoF, PostMdl.DoF);
end

% Checking Means
M1=mean(Sigmadraws,3);
M2=PostMdl.Omega/(PostMdl.DoF-nvar-1);
M3=PostMdl.Covariance;
for i=1:nd
    Sigmadraws(:,:,i)=iwishrnd(PostMdl.Omega, PostMdl.DoF);
end
M4=mean(Sigmadraws,3);

% Display results
fprintf('Distance from M1 to M3: %.4f\n', norm(M1 - M3, 'fro'));
fprintf('Distance from M2 to M3: %.4f\n', norm(M2 - M3, 'fro'));
fprintf('Distance from M4 to M3: %.4f\n', norm(M4 - M3, 'fro'));

%% Model definition using bayesvarm + ModelType ='conjugate'

disp('bayesvarm + ModelType = conjugate')

PriorMdl_bayesvarmconjugate = bayesvarm(nvar,nlag,ModelType='conjugate'); 

PriorMdl_bayesvarmconjugate.Omega=PphiBar;
PriorMdl_bayesvarmconjugate.DoF=nnuBar;
PriorMdl_bayesvarmconjugate.Mu=PpsiBar;
PriorMdl_bayesvarmconjugate.V= OomegaBar;
IncludeConstant=nex;                                
IncludeTrend=false; 

PostMdl_bayesvarmconjugate = PriorMdl_bayesvarmconjugate.estimate(Y, Display = 'off');

[~,Sigmadraws_bayesvarmconjugate] = simulate(PriorMdl_bayesvarmconjugate, Y, NumDraws = nd);

PostMdl_bayesvarmconjugate.Omega;
PostMdl_bayesvarmconjugate.DoF;

if T + PriorMdl_bayesvarmconjugate.DoF == PostMdl_bayesvarmconjugate.DoF
    fprintf('Degrees of freedom match: %d + %d = %d\n', T, PriorMdl_bayesvarmconjugate.DoF, PostMdl_bayesvarmconjugate.DoF);
else
    error('Degrees of freedom mismatch: %d + %d ≠ %d', T, PriorMdl_bayesvarmconjugate.DoF, PostMdl_bayesvarmconjugate.DoF);
end

% Checking Means
M5=mean(Sigmadraws_bayesvarmconjugate,3);
M6=PostMdl_bayesvarmconjugate.Omega/(PostMdl_bayesvarmconjugate.DoF-nvar-1);
M7=PostMdl_bayesvarmconjugate.Covariance;
for i=1:nd
    Sigmadraws(:,:,i)=iwishrnd(PostMdl_bayesvarmconjugate.Omega, PostMdl_bayesvarmconjugate.DoF);
end
M8=mean(Sigmadraws,3);

% Display results
fprintf('Distance from M5 to M3: %.4f\n', norm(M5 - M3, 'fro'));
fprintf('Distance from M6 to M3: %.4f\n', norm(M6 - M3, 'fro'));
fprintf('Distance from M7 to M3: %.4f\n', norm(M7 - M3, 'fro'));
fprintf('Distance from M8 to M3: %.4f\n', norm(M8 - M3, 'fro'));

PostMdl.Omega







