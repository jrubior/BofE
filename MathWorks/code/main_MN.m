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
nlag      =1;               % number of lags
nex       = 1;               % set equal to 1 if a constant is included; 0 otherwise
m         = nvar*nlag + nex; % number of exogenous variables
T         = size(Y,1) - nlag; % sample
pos       = []; % if there is variables that are stationary.

%% MN Prior Parameters

disp('MN Prior')

Vc=1e4;
lambda1=0.2;
lambda3=4;

% nnuBar
nnuBar = nvar+2;

% PphiBar
% Residual variance of AR(1) for each variable
diagPphiBar=zeros(nvar,1);
for i=1:nvar
    Mdl = arima(nlag, 0, 0);       % AR(1) model: p=1, d=0, q=0
    EstMdl = estimate(Mdl, Y(:,i),'Display','off');
    diagPphiBar(i)=EstMdl.Variance*(nnuBar-nvar-1);
end
PphiBar=diag(diagPphiBar);

% PpsiBar
PpsiBar= zeros(m,nvar);
diagPpsiBar=ones(nvar,1);
diagPpsiBar(pos)=0; % if there is variables that are stationary.
PpsiBar(1:nvar,:)=diag(diagPpsiBar);

% OomegaBar
OomegaBar = zeros(m, m);
if nex==1
    OomegaBar(end,end) = Vc;
end

KK=((nnuBar-nvar-1)*(lambda1^2))*inv(diag(((1:nlag)'.^lambda3)));
OomegaTilde=kron(KK,PphiBar\eye(nvar));
OomegaBar(1:end-nex, 1:end-nex) = OomegaTilde;

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





