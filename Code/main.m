%==========================================================================
%% housekeeping
%==========================================================================

clear variables;close all;userpath('clear');restoredefaultpath;clc;
tic;
rng('default'); % reinitialize the random number generator to its startup configuration
seed = 0; % original
rng(seed,'twister');  % set seed
currdir=pwd;

% Rutas de funciones auxiliares
addpath([currdir,'/helpfunctions']); 
addpath([currdir,'/init_parameters']); 
addpath([currdir,'/helpfunctions/ChrisSimsOptimize']);
addpath([currdir,'/helpfunctions/subroutines']); 
addpath([currdir,'/helpfunctions/ess_within_smc']); 

% Agregamos la ruta de datos y todas sus subcarpetas
cd(currdir)

%==========================================================================
%% load the Alpha data
%==========================================================================
load('data/data.mat')

%==========================================================================
%% choice of estimation sample, constant or varying volatility, and forecasting period
%==========================================================================

% --- Define estimation sample 
DATE_START = '1990-03-01'; 
DATE_END   = '2025-12-01';
DATE_END_DATA  = '2025-12-01';

Dates = datetime(Dates, 'InputFormat', 'dd-MMM-yyyy');

T0 = find(Time == datenum(DATE_START));
T1estim = find(Time == datenum(DATE_END));
Tend = find(Time == datenum(DATE_END_DATA));

% Covid Reference
T2019Q4 = find(Time == datenum(2019,12,01));
Tcovid = T2019Q4 - T0 + 2; % Se usará si hyperparameters = 'LP'
covid_periods = 3; % Se usará si hyperparameters = 'LP'

VarSelection = [1:11,13:14,16:17];

date = datetime(Dates(T1estim));

num = y(T0:T1estim,VarSelection);
VarNames = VarNames(VarSelection);

%=========================================================================
%% model setup
%==========================================================================
pos = find(VarNames == "VIX");    % stationary variables
nlag      = 5;                     % number of lags
nvar      = size(num,2);           % number of endogenous variables
nex       = 1;                     % 12 because of seasonal dummies
m         = nvar*nlag + nex;       % number of exogenous variables
horizon   = 35;                    % maximum horizon for IRFs
horizons  = [0,1,2,3,4,inf];       % horizons upon which sign and zero restrictions can be imposed
NS        = 1 + numel(horizons);   % number of objects in F(THETA) to which we impose sign and zero restrictions: F(THETA)=[A_0;L_{0};L_{1};L_{2};L_{3};L_{4};L_{inf}]
e         = eye(nvar);             % create identity matrix

M0  = 510000;                      % total iterations
Nburn  = 10000;                    % iterations of burn
save_every = 10;                   % save every save_every
nsave      = (M0-Nburn)/save_every; % we store nsave elements

conjugate  = 'none';
iter_show  = 1e3;
fixed_rf   = 0;
prior_only = 0;
label_R = 'cmy';
prior_type = 'minnesota';
hyperparameters = 'LP';
use_soc = 1;            % 1 = activar SOC, 0 = desactivar SOC (sum-of-coefficients)
use_dio = 1;            % 1 = activar DIO, 0 = desactivar DIO (dummy of initial observations)

if fixed_rf==1
    draw_Sigma=0;
    draw_B=0;
else
    draw_Sigma=1;
    draw_B=1;
end

elliptical_Q = 1; % 0 when using Chan, Matthes, Yu, 1 when using elliptical

%============================================================
%% Initialization setup
%==========================================

initialize = 'normal'; % 'normal' or 'none'
init_params_name = 'init_params_199003-202512_soc1_dio1_lag5_3shocks.mat' ; % provide an inizialization's name if initialize == 'none' 

myCluster = parcluster('Processes');
delete(myCluster.Jobs)

%=========================================================================
%% Tcovid based on hyperparameters estimation methodology
%=========================================================================

if strcmp(hyperparameters, 'GLP')
    Tcovid = [];
    covid_periods = 0; 
end

% Check de sanidad: Si usamos LP pero la muestra es PRE-COVID
if strcmp(hyperparameters, 'LP') && T1estim <= T2019Q4
    warning('hyperparameters="LP" pero la muestra termina antes de COVID. Forzando Tcovid = [].');
    Tcovid = [];
    covid_periods = 0;
end

%=========================================================================
%% Set up for organization of results 
%=========================================================================

% Quitar guiones de las fechas
dates_str = sprintf('%s-%s', ...
    datestr(datetime(DATE_START),'yyyymm'), ...
    datestr(datetime(DATE_END),'yyyymm'));
    
% Indicador SOC
if use_soc == 1
    soc_str = "soc1";
else
    soc_str = "soc0";
end

% Indicador DIO
if use_dio == 1
    dio_str = "dio1";
else
    dio_str = "dio0";
end
    
% Número de rezagos
lag_str = "lag" + string(nlag);

% Construir nombre final
sample = dates_str + "_" + soc_str + "_" + dio_str + "_" + lag_str;

%==========================================================================
%% Setup info/settings and IDENTIFYING RESTRICTIONS
%==========================================================================
info = SetupInfo(nvar,m, nlag,horizons,@(x)chol(x));
info.pos = pos;
info.use_soc = use_soc;
info.use_dio = use_dio;

%==========================================================================
%% Mappings
%==========================================================================
fo                 = @(x)f_h(x,info);
fo_inv             = @(x)f_h_inv(x,info);
fo_str2irfs        = @(x)StructuralToIRF(x,info);
fo_str2irfs_inv    = @(x)IRFToStructural(x,info);
gs_qr    = @(x)qr_unique(x);

%==========================================================================
%% write data in Rubio, Waggoner, and Zha (RES 2010)'s notation
%==========================================================================
y0bar = mean(num(1:nlag,:),1); % useful for GLP prior
info.y0bar = y0bar;

Tnum = size(num,1);
yt   = num(nlag+1:end,:);
T    = size(yt,1);
Tcovid = Tcovid-nlag;
xt = zeros(T,nvar*nlag+nex);

for i=1:nlag
    xt(:,nex+nvar*(i-1)+1:nex+nvar*i) = num((nlag-(i-1)):end-i,:);
end

if nex==1
    xt(:,1) = ones(T,1);
end

% write data in Zellner (1971, pp 224-227) notation
Y = yt; % T by nvar matrix of observations
X = xt; % T by (nvar*nlag+1) matrix of regressors

Ynd = Y; %% we will use Ynd for transformations, and keep Y "pure"
Xnd = X;

%==========================================================================
%% Prior for reduced-form parameters
%==========================================================================

switch prior_type
    case 'flat'

    %% prior
    nnuBar             = 0;
    OomegaBarInverse   = zeros(m);
    mmuBar             = zeros(m,nvar);  % Psi in the paper
    PpsiBar            = zeros(nvar);    % Phi in the paper

    case 'minnesota'

        SS=zeros(nvar,1);
        for i=1:nvar
            Tend=T; 
            if ~isempty(Tcovid)
                Tend=Tcovid-1; 
            end
            yols=Y(2:Tend,i);
            xols=[ones(Tend-1,1),Y(1:Tend-1,i)];
            bbetaols = (xols'*xols)\(xols'*yols);
            eols = yols-xols*bbetaols;
            SS(i)=eols'*eols/(size(yols,1)-size(xols,2));
        end
        MIN.ppsi  = SS./100;
        MAX.ppsi  = SS.*100;
        ppsi      = SS; %

        llambda0 = 0.2;     % Default
        alpha0   = 2;       % Default
        ppsi0    = SS;      % Default (LP usa esto, GLP lo sobrescribe)
        eta0     = [];      % Default (GLP usa esto, LP lo sobrescribe)
        info.Vc  = 10^7;    % Varianza del prior

        switch hyperparameters

            case 'fixed'

                llambda0 = 0.076676462;
                ppsi0    = ppsi; %zeros(nvar,1);
                alpha0 = 0.824764951;
                miu0 = 1;     
                theta0 = 1; 
            
            case 'GLP'

                % Giannone, Lenza, Primiceri (2015) 
                mode.llambda = .2;
                sd.llambda = .4;
                scalePSI   = 0.02^2; 
                mode.miu = 1;
                sd.miu =1; 
                mode.theta = 1;
                sd.theta =1;

                % transform mode and std into Gamma(k,theta) coefficients
                info.priorcoef.llambda= GammaCoef(mode.llambda,sd.llambda,0);
                info.priorcoef.alpha.PSI=scalePSI;
                info.priorcoef.beta.PSI=scalePSI;
                info.priorcoef.miu= GammaCoef(mode.miu,sd.miu,0);
                info.priorcoef.theta=GammaCoef(mode.theta,sd.theta,0);

                % bounds for maximization
                MIN.llambda = 0.0001;
                MAX.llambda = 5;
                MIN.alpha = 0.1;
                MAX.alpha = 5;
                MIN.miu = 0.0001;
                MAX.miu = 50;
                MIN.theta = 0.0001;
                MAX.theta = 50;

                % initial guess for hyper-parameters
                llambda  = 0.5;
                alpha = 2;
                miu = 1; 
                theta = 1; 
                
                inllambda = -log((MAX.llambda-llambda)./(llambda-MIN.llambda));
                inppsi    = -log((MAX.ppsi-ppsi)./(ppsi-MIN.ppsi));
                inalpha = -log((MAX.alpha-alpha)./(alpha-MIN.alpha));
                x0 = [inllambda;inppsi;inalpha];
                
                if use_soc == 1 
                    inmiu = -log((MAX.miu-miu)./(miu-MIN.miu));
                    x0 = [x0; inmiu];
                end 

                if use_dio == 1 
                    intheta = -log((MAX.theta-theta)./(theta-MIN.theta));
                    x0 = [x0; intheta];
                end 

                crit                              = 1e-16;
                nit                               = 1000;
                H0                                = eye(size(x0,1))*10;

                [fh,xh,~,~,itct,~,~]              = csminwel('negative_loglike_obj_mp_with_alpha',x0,H0,[],crit,nit,info,Xnd,Ynd,MIN,MAX,SS, covid_periods, Tcovid);

                % --------- parámetros base -----------------------------------
                idx = 1;
                
                llambda0 = MIN.llambda + (MAX.llambda-MIN.llambda)./ (1 + exp(-xh(idx, 1)));
                idx = idx + 1;
                
                ppsi0 = MIN.ppsi + (MAX.ppsi-MIN.ppsi)./ (1 + exp(-xh(idx:idx+info.nvar-1,1)));
                idx = idx + info.nvar;
                
                alpha0 = MIN.alpha + (MAX.alpha-MIN.alpha)./ (1 + exp(-xh(idx, 1)));
                idx = idx + 1;
                
                % --------- SOC ------------------------------------------------
                if use_soc == 1
                    miu0 = MIN.miu + (MAX.miu-MIN.miu)./ (1 + exp(-xh(idx,1)));
                    idx = idx + 1;
                end
                
                % --------- DIO ------------------------------------------------
                if use_dio == 1
                    theta0 = MIN.theta + (MAX.theta-MIN.theta)./ (1 + exp(-xh(idx,1)));
                    idx = idx + 1;
                end

                
                if itct > 1 % we re-set the random seed after using csminwel to guarantee the
                % posterior draws are identical across operating systems.
                    rng(seed,'twister');
                end

            case 'LP'
                mode.llambda = .2;
                sd.llambda = .4;
                mode.miu= 1;
                sd.miu=1;
                mode.theta= 1;
                sd.theta=1;
                info.priorcoef.llambda = GammaCoef(mode.llambda,sd.llambda,0);
                info.priorcoef.miu = GammaCoef(mode.miu,sd.miu,0);
                info.priorcoef.theta = GammaCoef(mode.theta,sd.theta,0);
                mode.eta(4) = .8;
                sd.eta(4) = .2;
                mosd = [mode.eta(4) sd.eta(4)];
                fun = @(x) BetaCoef(x,mosd);
                albet = fsolve(fun,[2,2])';
                info.priorcoef.eta4.alpha=albet(1);    
                info.priorcoef.eta4.beta=albet(2);

                % Bounds
                MIN.llambda = 0.0001; 
                MAX.llambda = 5;
                MIN.alpha = 0.1;   
                MAX.alpha = 5;
                MIN.miu = 0.0001;
                MAX.miu=50;
                MIN.theta = 0.0001;
                MAX.theta=50;
                
                % Guesses iniciales
                llambda  = 0.2;
                alpha = 2;
                miu= 1;
                theta = 1; 
                
                inllambda = -log((MAX.llambda-llambda)./(llambda-MIN.llambda));
                inalpha = -log((MAX.alpha-alpha)./(alpha-MIN.alpha));
                
                if isempty(Tcovid)

                    error('Running LP but Tcovid is empty. Cannot continue.')
                else
                    % --- CASO: LP con fecha COVID ---
                    
                    aux = mean(abs(Y(Tcovid:max([Tcovid+1 T]),:)-Y(Tcovid-1:max([Tcovid+1 T])-1,:))',1)./...
                          mean(mean(abs(Y(2:Tcovid-1,:)-Y(1:Tcovid-2,:))));
                    
                    if isempty(aux) || numel(aux) < covid_periods
                        eta0_initial = [ones(covid_periods,1); .8]; % Fallback guess
                    else
                        eta0_initial = [aux(1:covid_periods)'; .8]; % Guess inicial
                    end

                    % Escalas por periodo COVID
                    MIN.eta(1:covid_periods) = 1;
                    MAX.eta(1:covid_periods) = 500;

                    % Rho (decay)
                    MIN.eta(covid_periods+1) = 0.005;
                    MAX.eta(covid_periods+1) = 0.995;
                    
                    ineta = -log((MAX.eta'-eta0_initial)./(eta0_initial-MIN.eta'));

                    x0 = [inllambda;ineta;inalpha];
                    
                    if use_soc == 1 
                        inmiu = -log((MAX.miu-miu)./(miu-MIN.miu));
                        x0 = [x0; inmiu];
                    end 
    
                    if use_dio == 1 
                        intheta = -log((MAX.theta-theta)./(theta-MIN.theta));
                        x0 = [x0; intheta];
                    end 

                    % f0 = negative_loglike_obj_mp_with_alpha(x0, info, Xnd, Ynd, MIN, MAX, SS, covid_periods, Tcovid);
                    
                    crit = 1e-16; 
                    nit = 1000; 
                    H0 = eye(size(x0,1))*10;

                    [fh,xh,~,~,itct,~,~] = csminwel('negative_loglike_obj_mp_with_alpha',x0,H0,[],crit,nit,info,Xnd,Ynd,MIN,MAX, SS, covid_periods, Tcovid);
               
                    % número de parámetros eta
                    ncp = covid_periods + 1;
                    
                    % --------- parámetros base -----------------------------------
                    idx = 1;
                    
                    % llambda (primer elemento)
                    llambda0 = MIN.llambda + (MAX.llambda - MIN.llambda) ./ (1 + exp(-xh(idx,1)));
                    idx = idx + 1;
                    
                    % eta: ncp elementos (sustituye a ppsi)
                    eta_idx = idx : (idx + ncp - 1);
                    eta0 = MIN.eta(:) + (MAX.eta(:) - MIN.eta(:))./ (1 + exp(-xh(eta_idx)));
                    idx = idx + ncp;

                    % alpha (siguiente elemento)
                    alpha0 = MIN.alpha + (MAX.alpha - MIN.alpha) ./ (1 + exp(-xh(idx,1)));
                    idx = idx + 1;
                    
                    % --------- SOC ------------------------------------------------
                    if exist('use_soc','var') && use_soc == 1
                        miu0 = MIN.miu + (MAX.miu - MIN.miu) ./ (1 + exp(-xh(idx,1)));
                        idx = idx + 1;
                    end
                    
                    % --------- DIO ------------------------------------------------
                    if exist('use_dio','var') && use_dio == 1
                        theta0 = MIN.theta + (MAX.theta - MIN.theta) ./ (1 + exp(-xh(idx,1)));
                        idx = idx + 1;
                    end
                    
                    if itct > 1, rng(seed,'twister'); end
                end
        end 
 
        % --- 4. Re-ponderar Datos (Solo si hay COVID) ---
        if ~isempty(Tcovid) && ~isempty(eta0) && ~isempty(covid_periods)
        
            fprintf('--- Applying time-varying volatility weights---\n');
            invweights = ones(T,1);
     
            % Safety: truncar covid_periods si excede la muestra
            % disponible, o sea si pusimos 10 covid periods pero post covid
            % hay solo 5 periodos 
            covid_periods_eff = min(covid_periods, T - Tcovid + 1);
            if covid_periods_eff < covid_periods
                warning('covid_periods truncated to fit sample length.');
            end
        
            % 1) Bloque COVID
            covid_idx = Tcovid : (Tcovid + covid_periods_eff - 1);
       
            if numel(eta0) < covid_periods_eff + 1
                error('eta0 debe tener al menos covid_periods_eff + 1 elementos (último = rho).');
            end
        
            invweights(covid_idx) = eta0(1:covid_periods_eff);
        
            % 2) Decay post-COVID
            rho = eta0(end);
            post_start = Tcovid + covid_periods_eff;   % primer periodo después del bloque COVID
            if post_start <= T
                eta_last = eta0(covid_periods_eff);    % último multiplicador aplicado dentro del bloque COVID
                k = 1 : (T - post_start + 1);          % empieza en 1 para NO repetir eta_last
                invweights(post_start:T) = 1 + (eta_last - 1) * (rho .^ k)';
            end
        
            % Aplicar weights
            W = 1 ./ invweights;
            Ynd = diag(W) * Ynd;
            Xnd = diag(W) * Xnd;
        
        end
        Td=0;

        % ---- SOC ----
        if info.use_soc == 1 
            scl = 1/miu0;                                 % convención: (1/miu)
            ydsoc = scl * diag(info.y0bar);               % n x n diagonal
            if isfield(info,'pos') && ~isempty(info.pos)
                ydsoc(info.pos, info.pos) = 0;            % respetar pos si aplica
            end
            xdsoc = [zeros(info.nvar,1), scl*repmat(diag(info.y0bar),1,info.nlag)]; % n x (nex + nvar*nlag)
            Ynd = [Ynd; ydsoc];
            Xnd = [Xnd; xdsoc];
            Td = Td + size(Ynd,2);  
        end

        % ---- DIO ----
        if info.use_dio == 1
    
            scl = 1/theta0;   
            yddio = (scl * info.y0bar(:))';   % forzamos fila
            if isfield(info,'pos') && ~isempty(info.pos)
                yddio(info.pos) = 0;            % respetar 'pos' si aplica
            end
            xddio = [scl, scl * repmat(info.y0bar(:)', 1, info.nlag)];

            % anexar exactamente como hace sur
            Ynd = [Ynd; yddio];
            Xnd = [Xnd; xddio];

            % comportamiento idéntico al if sur...: asignar Td = 1
            Td = Td + 1;

       end 

        % --- 5. Configuración Final del Prior  ---
        nnuBar              = info.nvar + 2;
        PpsiBar             = diag(ppsi0); 
        oomegaBar=zeros(info.m,1);
        oomegaBar(1:info.nex,1)=info.Vc; 
        for ell=1:nlag
            oomegaBar(info.nex+1+(ell-1)*info.nvar:info.nex+ell*info.nvar,1) =  (llambda0^2)*(nnuBar-info.nvar-1)./((ell^alpha0)*ppsi0);
        end
        OomegaBar=diag(oomegaBar);
        OomegaBarInverse=diag(1./oomegaBar);
        mmuBar = zeros(info.m,info.nvar);
        diagmmuBar=ones(info.nvar,1);
        diagmmuBar(info.pos)=0;   
    
        mmuBar(info.nex+1:info.nvar+info.nex,:)=diag(diagmmuBar);   
end

info.lambda = llambda0;
info.psi = ppsi0;
info.alpha = alpha0;
info.eta = eta0;

if use_soc == 1
    info.miu = miu0;
end

if use_dio == 1
    info.dio = theta0;
end

info.nnuBar = nnuBar;
info.PpsiBar = PpsiBar;
info.oomegaBar = OomegaBar;
info.oomegaBarInverse = OomegaBarInverse;
info.mmuBar = mmuBar;
info.diagmmuBar = diagmmuBar;

%==========================================================================
%% Posterior for reduced-form parameters
%==========================================================================

nnuTilde            = T + nnuBar + Td;
OomegaTilde         = (Xnd'*Xnd  + OomegaBarInverse)\eye(m);
OomegaTildeInverse  =  Xnd'*Xnd  + OomegaBarInverse;
mmuTilde            = (Xnd'*Xnd  + OomegaBarInverse)\(Xnd'*Ynd + OomegaBarInverse*mmuBar);
PpsiTilde           = Ynd'*Ynd + PpsiBar + mmuBar'*OomegaBarInverse*mmuBar - mmuTilde'*OomegaTildeInverse*mmuTilde;
PpsiTilde           = (PpsiTilde'+PpsiTilde)*0.5;

info.nnuTilde =nnuTilde;
info.OomegaTilde=OomegaTilde;
info.OomegaTildeInverse=OomegaTildeInverse;
info.mmuTilde =mmuTilde ;
info.PpsiTilde=PpsiTilde;
info.pred=info.m;

%==========================================================================
%% Parameters to fix reduced-form
%==========================================================================

BHat = mmuTilde; %(X'*X)\(X'*Y); %
SigmaHat =  PpsiTilde/(nnuTilde+size(PpsiTilde,1)+1); %(Y-X*BHat)'*(Y-X*BHat)/size(Y,1); %
SigmaHat = (SigmaHat'+SigmaHat)*0.5;

%==========================================================================
%% Load restrictions from excel
%==========================================================================

%IMPORTANT: if initialize == none, it is important that excel restrictions
%match restrincitions from initial parameters 

excelPath = fullfile(currdir, 'restrictions.xlsx');

% IRFs restrictions
IRF_sheetSigns = 'IRF_signs';
IRF_sheetRanks = strings(0); % 'IRF_rank' or strings(0)

% A0 restrictions
A0_sheetSigns  = 'A0_signs';   %      coef > a  
A0_sheetBounds = 'A0_bounds';  %      valor de 'a'

info.IRFepsilon = 0.01; %  -epsilon <   impacto   < epsilon
info.A0epsilon  = 0.01; %  -epsilon < coeficiente < epsilon

% Load restrictions as info.Ss
info = read_restrictions(info, excelPath, IRF_sheetSigns, ...
                                          IRF_sheetRanks, ...
                                          A0_sheetSigns, ...
                                          A0_sheetBounds);


%==========================================================================
%% Starts initiliazation
%==========================================================================
init = 0;
Sigmadraw=SigmaHat;
Bdraw=BHat;
  
% =====================================
%%  New initialization
% =====================================
switch initialize

    case 'normal'

        disp('initilization of Q ...');
        
        init_params.Q_init = [];

        init_time =tic;

        [Q0,~]      = qr(randn(info.nvar));
        Svec = SRvec(Bdraw, Sigmadraw, Q0, info);

        tStart = tic;
        TIME_INIT = toc(tStart);
        function_restrictions_Svec = @SRvec;

        rst = ess_within_smc_init(Ynd,function_restrictions_Svec,info, init_params);
        TIME_INIT = toc(tStart);

        z_old_B = vec(squeeze(rst.mat_B(:,:,1)));
        z_old_R = rst.mat_R(:,1);
        R_init = z_old_R;
        z_old_X = rst.mat_X(:,1);
        X_init = z_old_X;
        Bdraw    = rst.mat_B(:,:,1);
        B_init = Bdraw;
        Sigmadraw = rst.mat_S(:,:,1);
        S_init = Sigmadraw;
        Qdraw = rst.mat_Q(:,:,1);
        Q_init = Qdraw;
        disp('initilization of Q ... done ...');

        time_to_initialize = toc(init_time);
        disp(time_to_initialize)
        init_params_save = fullfile('init_parameters',"init_params_" + sample + ".mat");
        save(init_params_save, 'z_old_B', 'z_old_X', 'z_old_R', 'B_init', 'S_init', 'Q_init', 'R_init', 'X_init');

    case 'none'

        disp('Using Q and params from old initialization directly');

        init_params = load(init_params_name);
        
        z_old_R = init_params.R_init;
        R_init = z_old_R;

        z_old_X = init_params.X_init;
        X_init = z_old_X;

        Bdraw    = init_params.B_init;
        B_init = Bdraw;
        z_old_B = vec(squeeze(Bdraw));

        Sigmadraw = init_params.S_init;
        S_init = Sigmadraw;

        Qdraw = init_params.Q_init;
        Q_init = Qdraw;

        disp('initilization of Q ... done ...');

end

%==========================
%% Store draws before Gibbs 
%==========================
% definitions used to store orthogonal-reduced-form draws, volume elements, and unnormalized weights

Bdraws         = cell([nsave,1]); % reduced-form lag parameters
Sigmadraws     = cell([nsave,1]); % reduced-form covariance matrices
Qdraws         = cell([nsave,1]); % orthogonal matrices
A0draws        = cell([nsave,1]); % structural contemporaneous parameters

%==========================
%% Definitions to facilitate the draws from B|Sigma
%==========================
hh              = info.h;
cholOomegaTilde = hh(OomegaTilde)'; % this matrix is used to draw B|Sigma below

chol_lower_OomegaTilde = chol(OomegaTilde, 'lower');
chol_lower_OomegaTilde_inv = eye(size(OomegaTilde))/chol_lower_OomegaTilde;

if prior_only==1
    cholOomegaBar   = hh(OomegaBar)'; % this matrix is used to draw B|Sigma below
end

Rmean = nan(nvar,nnuTilde);
for i=1:nnuTilde
    Rmean(:,i) = zeros(nvar,1);
end
Rvariance = zeros(nvar*nnuTilde,nvar*nnuTilde);
for i=1:nnuTilde
    Rvariance((i-1)*nvar+1:i*nvar,(i-1)*nvar+1:i*nvar) = inv(PpsiTilde);
end

chol_lower_Rvariance = chol(Rvariance,'lower');
chol_lower_Rvariance_common = chol(inv(PpsiTilde), 'lower');		

SigmaHatInv = inv(SigmaHat);
R_base = chol(SigmaHatInv, 'lower');

%==========================================================================
%% Initialize counters to track state of computations 
%==========================================================================

counter = 1;
count   = 0;
record = 1;
record_save = 0;
save_indicator = 0;

switch label_R

    case 'cmy'

        function_restrictions = @restriction_baseline;
end

% check restrictions

function_restrictions_i = @(ag0, ag1, ag2, ag3) function_restrictions(ag0, ag1, ag2, ag3);
S_prop = function_restrictions_i(Bdraw, Sigmadraw, Qdraw, info);

if S_prop == 0

    disp("WARNING!!!! Restrictions not checked - PROBLEM WITH INITIALIZATION -")
end

tstart = tic;

%======================================================================
%% Starts GIBBS SAMPLER 
%======================================================================

ind_save=1; 
hh=@(x)chol(x);
info.hh=hh;
model=info;
model.save_every=save_every; 

while record<=M0
    
    % ===================================================================
    % BLOQUE: SLICE SAMPLERS
    % ===================================================================

    %% setting function restrictions
    function_restrictions_i = @(ag0, ag1, ag2, ag3) function_restrictions(ag0, ag1, ag2, ag3);

    %% Draw Q
    slice.scale_z    = 1;
    slice.mean       = zeros(model.nvar*model.nvar,1);
    slice.chol_cov_z = chol(eye(model.nvar*model.nvar),'lower');
    slice.chol_cov_z = 1; 
    slice.fcn_lik    = @(z_prop) loglike_Q(z_prop,function_restrictions_i,gs_qr,fo_inv,fo_str2irfs,Bdraw,Sigmadraw,info);
    slice.nobs       = model.nvar*model.nvar;
    lik_old = slice.fcn_lik(z_old_X);
    [z_old_X, ~, n_try] = slice_sampling_v02(slice, z_old_X, lik_old);
    [Qdraw,Rdraw] =  gs_qr(reshape(z_old_X,model.nvar,model.nvar));

    %% Draw Sigma
    slice.scale_z    = 1;
    slice.mean       = vec(Rmean);
    slice.chol_cov_z_common = chol_lower_Rvariance_common;
    slice.nobs1 = model.nvar;
    slice.nobs2 = nnuTilde;
    slice.fcn_lik    = @(z_prop) loglike_Sigma_mc_v02(z_prop,function_restrictions_i,fo_inv,fo_str2irfs,Bdraw,Qdraw,model,model.nnuTilde,model.mmuTilde,chol_lower_OomegaTilde,chol_lower_OomegaTilde_inv);
    slice.nobs       = model.nvar*nnuTilde;
    lik_old = slice.fcn_lik(z_old_R);
    [z_old_R, ~, n_try] = slice_sampling_v02_Sigma(slice, z_old_R, lik_old);
    R_old = reshape(z_old_R,model.nvar,nnuTilde);
    Sigmadraw=inv(R_old*R_old');
   
    %% Draw B
    cholLSigmadraw = model.hh(Sigmadraw)';
    slice.scale_z    = 1;
    slice.mean       = vec(mmuTilde);    
    slice.chol_cov_z1 = cholLSigmadraw;%chol(Rvariance,'lower');
    slice.chol_cov_z2 = cholOomegaTilde;%chol(Rvariance,'lower');  
    slice.nobs1 = model.nvar;
    slice.nobs2 = model.m;
 
    slice.fcn_lik    = @(z_prop) loglike_B(z_prop,function_restrictions_i,fo_inv,fo_str2irfs,Sigmadraw,Qdraw,model,model.mmuTilde,model.OomegaTilde);
    slice.nobs       = model.nvar*model.m;
    lik_old = slice.fcn_lik(z_old_B);
  
    [z_old_B, ~, n_try] = slice_sampling_v02_kron(slice, z_old_B, lik_old);
    Bdraw = reshape(z_old_B,model.m,model.nvar);

    % ===================================================================
    % GUARDAR SORTEOS (DESPUÉS DEL BURN-IN)
    % ===================================================================
    % Mueve 'count=count+1' aquí para contar solo las iteraciones válidas
    count=count+1; 
    
    if record > Nburn && mod(record,save_every)==0 
        record_save = record_save + 1;
        save_indicator = 1;
    end

    if record > Nburn  && save_indicator == 1
        Bdraws{record_save,1}     = Bdraw;
        Sigmadraws{record_save,1} = Sigmadraw;
        Qdraws{record_save,1}     = Qdraw;
    end

    record=record+1;
    save_indicator = 0;

    if counter==iter_show
        display(['Number of draws completed = ',num2str(record-1)])
        display(['Remaining draws to store = ',num2str(M0-(record-1))])
        counter =0;
    end
    counter = counter + 1;

end

telapsed = toc(tstart);

%% Store draws: L, cumL and A0
L    = zeros(horizon+1,nvar,nvar,nsave);
cumL    = zeros(horizon+1,nvar,nvar,nsave);

for s=1:(nsave)

    Bdraw =     Bdraws{s,1} ;
    Sigmadraw = Sigmadraws{s,1} ;
    Qdraw=Qdraws{s,1};

    BSigmaQ = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];

     structpara = f_h_inv(BSigmaQ,info);
    n2                        = nvar*nvar;
    A0draws{s,1}    = reshape(structpara(1:n2),nvar,nvar);

    LIRF = IRF_horizons(structpara, nvar, nlag, m, nex, 0:horizon);

    for h=0:horizon
        L(h+1,:,:,s) =  LIRF(1+h*nvar:(h+1)*nvar,:);
      
        for i=1:nvar
            cumL(h+1,1:4,i,s)   = sum(L(1:h+1,1:4,i,s),1);
        end

    end

end

if exist('A0draws','var'),           A0draws = ensure3D(A0draws); end
if exist('Bdraws','var'),            Bdraws = ensure3D(Bdraws); end
if exist('Sigmadraws','var'),    Sigmadraws = ensure3D(Sigmadraws); end
if exist('Qdraws','var'),            Qdraws = ensure3D(Qdraws); end

% ***********
% Histograma
% ***********
normalizationVariables = ["USRR","USCPI","USGDP"];
outdir = fullfile(pwd,'/results/structural coefficients/');
[tf, denominators] = ismember(normalizationVariables, VarNames);
nShocks = size(normalizationVariables,2);
for col = 1:nShocks
    fig = figure;
    for i = 1:nvar
        
        coef = squeeze(A0draws(i,col,:) ./ -A0draws(denominators(col),col,:));
        coef = coef(isfinite(coef));
    
        subplot(4,4,i);
    
        if i ~= denominators(col)
            xl = prctile(coef,[2.5 97.5]);
        else
            xl = [0.90 1.1];
        end
    
        coef_plot = coef(coef >= xl(1) & coef <= xl(2));
    
        histogram(coef_plot, ...
            'Normalization','pdf', ...
            'NumBins',30);

        title(VarNames{i}, 'Interpreter','none');
        grid on;
    
    end
    sgtitle(['A0 col ' num2str(col) ' coeficients (normalized)']);
    filename = fullfile(outdir,sprintf('A0_col_%d.jpg',col));
    exportgraphics(fig,filename,'Resolution',300);
end

% ***********
% Tabla
% ***********
filename = fullfile(pwd,'/results/structural coefficients/A0_summary.xlsx');

for col = 1:nShocks

    sheetname = sprintf('Column_%d', col);

    results = [];

    for r = 1:nvar

        coef = squeeze(A0draws(r,col,:) ./ -A0draws(denominators(col),col,:));
        coef = coef(isfinite(coef));

        mediana = median(coef);
        p16     = prctile(coef,16);
        p25     = prctile(coef,25);
        p75     = prctile(coef,75);
        p84     = prctile(coef,84);
        minimo  = min(coef);
        maximo  = max(coef);

        results = [results;
                   r mediana minimo p16 p25 p75 p84 maximo];

    end

    T = array2table(results, ...
        'VariableNames', ...
        {'RowIdx','Median','Min','P16','P25','P75','P84','Max'});

    T.Variable = VarNames(T.RowIdx)';
    T = movevars(T,'Variable','Before','Median');

    writetable(T, filename, ...
        'Sheet', sheetname, ...
        'Range', 'A1');

end


%=========================================================================
%% Store results - extract quantiles and median -
% %==========================================================================

store_results = true;
if store_results == true

    matfiles_path = fullfile(currdir, 'results', 'matfiles');
    if ~exist(matfiles_path, 'dir')
        mkdir(matfiles_path); 
        disp(['Carpeta creada: ', matfiles_path]);
    end

    excel_path = fullfile(currdir, 'results');
    if ~exist(excel_path, 'dir')
        mkdir(excel_path);
    end

    % Definir nombres de archivo usando rutas completas

    % 1. Convertir el número de draws a string (para evitar el error de tipo)
    num_draws_str = string(num2str(int64((M0-Nburn)/1000)));

    % 2. Concatenar todas las partes directamente
    filename_base = "results_" + string(hyperparameters) + ...
                    "_estim_" + string(dates_str) + ...
                    "_" + num_draws_str + "k_Draws";
    % Construir nombre de archivo Excel usando el operador '+' para strings
    filename_excel_str = filename_base + ".xlsx"; 
    filename_excel = fullfile(excel_path, char(filename_excel_str)); 

percentiles = 10:10:90;
nVar = size(num,2);

%==========================================================================
%% Save Beta and Sigma into excel sheets and matlab matrix
%==========================================================================


    % Save Beta draws 
    matfile_name_beta = sample + "_beta_draws.mat"; 
    full_path_beta = fullfile(matfiles_path, char(matfile_name_beta));
    % save(full_path_beta,'Bdraws')

    mediana_params = median(Bdraws, 3);
    std_params     = std(Bdraws, 0, 3);
    combined = [mediana_params, nan(size(Bdraws,1),1), std_params];
    colNames = [ strcat("mediana_eq", string(1:nVar)), " ", strcat("sd_eq", string(1:nVar)) ];
    T = array2table(combined, 'VariableNames', colNames);
    writetable(T, filename_excel, 'Sheet', 'Beta');

    % Save Sigma draws
    matfile_name_sigma = sample + "_sigma_draws.mat";
    full_path_sigma = fullfile(matfiles_path, char(matfile_name_sigma));
    % save(full_path_sigma,'Sigmadraws') 

    mediana_params = median(Sigmadraws, 3);
    std_params     = std(Sigmadraws, 0, 3);
    combined = [mediana_params, nan(size(Sigmadraws,1),1), std_params];
    colNames = [ strcat("mediana_eq", string(1:nVar)), " ", strcat("sd_eq", string(1:nVar)) ];
    T = array2table(combined, 'VariableNames', colNames);
    writetable(T, filename_excel, 'Sheet', 'Sigma');

    % Hyperparameters
    p_lambda = info.lambda;
    p_alpha  = info.alpha;
    p_psi    = info.psi;
    
    % SOC
    if use_soc==1
        p_miu    = info.miu;
    else
        p_miu    = nan(1);
    end
    
    % DIO
    if use_dio==1
        p_theta    = info.dio;
    else
        p_theta    = nan(1);
    end

    if isempty(Tcovid)
        TNames = ["Lambda", "Alpha", "Psi"+string(1:nvar), "miu", "theta"];
        T = array2table([p_lambda', p_alpha', p_psi', p_miu', p_theta'], 'VariableNames', TNames);
    else
        p_eta =  info.eta;
        TNames = ['Lambda', 'Alpha', "Psi"+string(1:nvar), "miu", "theta", "Eta"+string(0:2), "Rho"];
        T = array2table([p_lambda', p_alpha', p_psi', p_miu', p_theta', p_eta'], 'VariableNames', TNames);
    end
    writetable(T, filename_excel, 'Sheet', 'Hiperparámetros');   

    % Save Time spent on running the algorithm
    TimeNames = ["Total time mins", "mins by 1000 iterations"];
    T = array2table([telapsed/60, 1000*telapsed/(M0*60)], 'VariableNames', TimeNames);
    writetable(T, filename_excel, 'Sheet', 'Time Elapsed'); 

    % Save IRFs y resultados generales
    savefile_name = sample + "_rgibbs_results.mat";
    full_path_savefile = fullfile(matfiles_path, char(savefile_name));
    save(full_path_savefile,'L','cumL','horizon','telapsed','-v7.3');

end

