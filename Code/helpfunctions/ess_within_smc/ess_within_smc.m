function rst = ess_within_smc(data_num, function_restrictions_Svec, varargin)
% ESS_WITHIN_SMC  SMC with Elliptical Slice Sampling under sign/zero restrictions
%
% SYNTAX
%   rst = ess_within_smc(data_num, function_restrictions_Svec)
%   rst = ess_within_smc(data_num, function_restrictions_Svec, Name, Value, ...)
%
% DESCRIPTION
%   Runs an SMC sampler for a VAR with (potentially) identification
%   restrictions, using elliptical slice sampling (ESS) as the mutation step.
%
% REQUIRED INPUTS
%   data_num : T-by-n numeric matrix of observables.
%   function_restrictions_Svec : function handle taking (B, Sigma, Q, info)
%       and returning a vector Svec of restriction values (>=0 if satisfied).
%
% NAME-VALUE ARGUMENTS (OPTIONS)
%   'nlag'          (default 4)  [nonnegative integer]
%       Number of VAR lags. Affects state dimension m = nvar*nlag + nex.
%       Larger values increase parameter dimension and computation.
%
%   'prior_type'    (default 'flat')  ['flat' | 'minnesota' | other strings will be supported soon]
%       Chooses the prior family. In this file, 'flat' implies:
%       nnuBar = 0, OomegaBarInverse = 0, mmuBar = 0, PpsiBar = 0.
%       Extend this to other priors as needed.
%
%   'smc_eff_zero'  (default 1e-2)  [positive scalar]
%       Effective “zero” tolerance for the restriction vector Svec in the
%       cutoff search. Smaller values enforce tighter satisfaction of
%       restrictions; too small may cause particle degeneracy.
%
%   'smc_target_ar' (default 0.25)  (0,1)
%       Target acceptance rate used by the adaptive cutoff search to set
%       the incremental tempering level. Smaller values make each SMC step
%       more aggressive (larger moves) but can reduce effective sample size.
%
%   'smc_nsim'      (default 10000)  [positive integer]
%       Number of particles. Improves Monte Carlo accuracy but increases
%       memory and runtime roughly linearly in smc_nsim.
%
%   'smc_max_iter'  (default 200)  [positive integer]
%       Maximum number of SMC recursion steps (temperatures). Acts as a
%       guardrail to stop the algorithm if the cutoff search hasn’t reached
%       the final target yet.
%
%   'smc_ntran'     (default 5)  [positive integer]
%       Number of ESS mutation transitions per SMC stage (per particle).
%       Higher values improve within-stage mixing at the cost of runtime.
%
%   'smc_early_stop' (default false) [logical]
%       If true, stop the SMC recursion early when we find draw(s) that
%       satisfy restrictions (useful for finding initial value)
% NOTES
%   * The algorithm adaptively searches for a cutoff (quantile) on Svec to
%     achieve the target acceptance rate. Particles are reweighted, resampled,
%     and mutated via ESS at each stage until the final target (smc_eff_zero)
%     is met or smc_max_iter is reached.
%   * Memory considerations: storing B/Sigmas/etc. for smc_nsim=10k can be
%     heavy; consider reducing smc_nsim for prototyping.
%
% EXAMPLES
%   % All defaults
%   rst = ess_within_smc(data_num, @restriction_Svec_arrw_2018_more_zeros);
%
%   % Custom VAR lag and more particles
%   rst = ess_within_smc(data_num, @restriction_Svec_arrw_2018_more_zeros, ...
%                        'nlag', 6, 'smc_nsim', 20000);
%
%   % Tighter restriction tolerance and more ESS transitions per stage
%   rst = ess_within_smc(data_num, @restriction_Svec_arrw_2018_more_zeros, ...
%                        'smc_eff_zero', 1e-3, 'smc_ntran', 10);


% hardcoded options (to be relaxed)
reoptimizeMH = 0;
seed = 1234;


%% ---------------- Options (with defaults) ----------------

p = inputParser;
p.FunctionName = 'ess_within_smc';

addParameter(p, 'nlag',          4,       @(x)isnumeric(x) && isscalar(x) && x>=0 && mod(x,1)==0);
addParameter(p, 'prior_type',    'flat',  @(s)ischar(s) || (isstring(s) && isscalar(s)));
addParameter(p, 'smc_eff_zero',  1e-2,    @(x)isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'smc_target_ar', 0.25,    @(x)isnumeric(x) && isscalar(x) && x>0 && x<1);
addParameter(p, 'smc_nsim',      10000,   @(x)isnumeric(x) && isscalar(x) && x>0 && mod(x,1)==0);
addParameter(p, 'smc_max_iter',  200,     @(x)isnumeric(x) && isscalar(x) && x>0 && mod(x,1)==0);
addParameter(p, 'smc_ntran',     5,       @(x)isnumeric(x) && isscalar(x) && x>0 && mod(x,1)==0);
addParameter(p, 'smc_early_stop', false,   @(x)islogical(x) && isscalar(x));
addParameter(p, 'info', []);
parse(p, varargin{:});
opt = p.Results;

%% Unpack options

% --- options for VAR
nlag = opt.nlag;
nvar = size(data_num,2);
nex  = 1;
m    = nvar*nlag + nex; % number of exogenous variables
prior_type = char(opt.prior_type);

% --- options for SMC
smc_eff_zero  = opt.smc_eff_zero;   % effective zero for restrictions
smc_target_ar = opt.smc_target_ar;  % target acceptance rate (smaller -> more aggressive)
smc_nsim      = opt.smc_nsim;       % number of particles
smc_max_iter  = opt.smc_max_iter;   % maximum number of iterations
smc_ntran     = opt.smc_ntran;      % number of transitions
smc_early_stop = opt.smc_early_stop; % early stop for just finding initial value

% --- options for ESS

%% Pack options into the structure
info = opt.info;
info.nvar = nvar;
info.npredetermined = m;
info.h = @(x)chol(x);
info.nlag = nlag;
info.nex = nex;
info.m = m;

%%  Write data in Rubio, Waggoner, and Zha (RES 2010)'s notation
%==========================================================================
% yt(t) A0 = xt(t) Aplus + constant + et(t) for t=1...,T;
% yt(t)    = xt(t) B     + ut(t)            for t=1...,T;
% x(t)     = [constant,yt(t-1), ... ,yt(t-nlag)];
% matrix notation yt = xt*B + ut;
% xt=[yt_{-1} ones(T,1)];
%==========================================================================
yt   = data_num(nlag+1:end,:);
T    = size(yt,1);

xt = zeros(T,nvar*nlag+nex);
for i=1:nlag
    xt(:,nex+nvar*(i-1)+1:nex+nvar*i) = data_num((nlag-(i-1)):end-i,:) ;
end

if nex==1
    xt(:,1)=ones(T,1);
end

% write data in Zellner (1971, pp 224-227) notation
Y = yt; % T by nvar matrix of observations
X = xt; % T by (nvar*nlag+1) matrix of regressors

Ynd = Y;
Xnd = X;


%% Useful objects to pre-compute

% --- prior
switch prior_type
    case 'flat'
        % --- flat prior
        nnuBar             = 0;
        OomegaBarInverse   = zeros(m);
        mmuBar             = zeros(m,nvar);  % Psi in the paper
        PpsiBar            = zeros(nvar);    % Phi in the paper
    case 'minnesota'

        % Giannone, Lenza, Primiceri (2015) hyperpriors modes
        mode.llambda = .2;

        % Giannone, Lenza, Primiceri (2015) hyperpriors std deviations
        sd.llambda = .4;

        % Giannone, Lenza, Primiceri (2015) scale and shape of the IG on psi/(d-n-1)
        scalePSI   = 0.02^2;


        % transform mode and std into Gamma(k,theta) coefficients
        info.priorcoef.llambda= GammaCoef(mode.llambda,sd.llambda,0);
        info.priorcoef.alpha.PSI=scalePSI;
        info.priorcoef.beta.PSI=scalePSI;
        info.pos = []; % stationary variables


        % bounds for maximization
        MIN.llambda = 0.0001;
        MAX.llambda = 5;



        % initial guess for hyper-parameters
        llambda  = 0.5;
        info.Vc  = 10^7;

        % initial guess for psi (phi in the paper)
        % residual variance of AR(1) for each variable
        SS=zeros(nvar,1);
        for i=1:nvar
            yols=Y(2:end,i);
            xols=[ones(T-1,1),Y(1:end-1,i)];
            bbetaols = (xols'*xols)\(xols'*yols);
            eols = yols-xols*bbetaols;
            SS(i)=eols'*eols/(size(yols,1)-size(xols,2));
        end
        MIN.ppsi  = SS./100;
        MAX.ppsi  = SS.*100;
        ppsi      = SS;
        inllambda = -log((MAX.llambda-llambda)./(llambda-MIN.llambda));
        inppsi    = -log((MAX.ppsi-ppsi)./(ppsi-MIN.ppsi));

        x0 = [inllambda;inppsi];


        f0 = negative_loglike_obj_mp(x0,info,Xnd,Ynd,MIN,MAX);

        crit                              = 1e-16;
        nit                               = 1000;
        H0                                = eye(size(x0,1))*10;


        if reoptimizeMH==1

            [fh,xh,~,~,itct,~,~]              = csminwel('negative_loglike_obj_mp',x0,H0,[],crit,nit,info,Xnd,Ynd,MIN,MAX);

            save('xh.mat','xh','fh');
        else
            load('xh.mat','xh','fh');
        end

        llambda0_xh = MIN.llambda+(MAX.llambda-MIN.llambda)./(1+exp(-xh(1,1)));
        ppsi0_xh    = MIN.ppsi+(MAX.ppsi-MIN.ppsi)./(1+exp(-xh(2:info.nvar+1,1)));%

        llambda0 =  llambda0_xh;
        ppsi0    =  ppsi0_xh;


        % we re-set the random seed after using csminwel to guarantee the
        % posterior draws are identical across operating systems.
        rng(seed,'twister');



        % set position where variables are in first differences

        % inverse Wishart prior parameters
        nnuBar              = info.nvar + 2;
        PpsiBar             = diag(ppsi0);
        oomegaBar=zeros(info.m,1);
        oomegaBar(1:info.nex,1)=info.Vc;
        for ell=1:info.nlag
            oomegaBar(info.nex+1+(ell-1)*info.nvar:info.nex+ell*info.nvar,1) =  (llambda0^2)*(nnuBar-info.nvar-1)./((ell^2)*ppsi0);
        end
        OomegaBar=diag(oomegaBar);
        OomegaBarInverse=diag(1./oomegaBar);
        mmuBar = zeros(info.m,info.nvar);
        diagmmuBar=ones(info.nvar,1);
        diagmmuBar(info.pos)=0;   % Set to zero the prior mean on the first own lag for variables selected in the vector pos

        mmuBar(info.nex+1:info.nvar+info.nex,:)=diag(diagmmuBar);


        Y = Ynd;
        X = Xnd;
        T = size(Y,1);


end

% --- preliminary calculation
nnuTilde            = T +nnuBar;
OomegaTilde         = (X'*X  + OomegaBarInverse)\eye(m);
OomegaTildeInverse  =  X'*X  + OomegaBarInverse;
mmuTilde            = OomegaTilde*(X'*Y + OomegaBarInverse*mmuBar);
PpsiTilde           = Y'*Y + PpsiBar + mmuBar'*OomegaBarInverse*mmuBar - mmuTilde'*OomegaTildeInverse*mmuTilde;
PpsiTilde           = (PpsiTilde'+PpsiTilde)*0.5;

% PpsiTilde = PpsiTilde + 10000 * eye(size(PpsiTilde));

cholPpsiTilde = chol(PpsiTilde,'lower');
inv_cholPpsiTilde_prime = inv(cholPpsiTilde');

hh              = info.h;
cholOomegaTilde = hh(OomegaTilde)'; % this matrix is used to draw B|Sigma below

% --- iid draw from the posterior (sitting here for testing and to get dimension)
[Bdraw, Sigmadraw, Rdraw, Qdraw, Xdraw] = iid_sampler(mmuTilde, cholOomegaTilde, inv_cholPpsiTilde_prime, nnuTilde, info);
Svec = function_restrictions_Svec(Bdraw, Sigmadraw, Qdraw, info);

% --- Setting up for ESS
% definition to facilitate the draws from B|Sigma
hh              = info.h;
cholOomegaTilde = hh(OomegaTilde)'; % this matrix is used to draw B|Sigma below

chol_lower_OomegaTilde = chol(OomegaTilde, 'lower');
chol_lower_OomegaTilde_inv = eye(size(OomegaTilde))/chol_lower_OomegaTilde;

% if prior_only==1
%     cholOomegaBar   = hh(OomegaBar)'; % this matrix is used to draw B|Sigma below
% end


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


% collect information for ess
info_ess = [];
info_ess.Rmean = Rmean;
info_ess.chol_lower_Rvariance_common = chol_lower_Rvariance_common;
info_ess.nnuTilde = nnuTilde;

info_ess.mmuTilde = mmuTilde;
info_ess.chol_lower_OomegaTilde = chol_lower_OomegaTilde;
info_ess.chol_lower_OomegaTilde_inv = chol_lower_OomegaTilde_inv;
info_ess.cholOomegaTilde = cholOomegaTilde;


%% Actual loop

% --- Initialization
mat_B = nan([size(Bdraw), smc_nsim]);
mat_S = nan([size(Sigmadraw), smc_nsim]);
mat_R = nan([numel(Rdraw), smc_nsim]);
mat_Q = nan([size(Qdraw), smc_nsim]);
mat_X = nan([numel(Xdraw), smc_nsim]);
mat_Svec = nan([numel(Svec), smc_nsim]);
mat_W = ones(1, smc_nsim); % initial weihgting is 1 because they came from prior

for i=1:1:smc_nsim

    % iid draws
    [Bdraw, Sigmadraw, Rdraw, Qdraw, Xdraw] = iid_sampler(mmuTilde, cholOomegaTilde, inv_cholPpsiTilde_prime, nnuTilde, info);
    Svec = function_restrictions_Svec(Bdraw, Sigmadraw, Qdraw, info);

    % store
    mat_B(:,:,i) = Bdraw;
    mat_S(:,:,i) = Sigmadraw;
    mat_R(:,i)   = Rdraw;
    mat_Q(:,:,i) = Qdraw;
    mat_X(:,i)   = Xdraw;
    mat_Svec(:,i) = Svec;

end

% --- Recursion
smc_iter = 2;
smc_continue = 1;

% *** while smc_continue
while smc_continue

    % --- (a) finding a good cut-off
    % logic: finding quantile of each Svec such that it gives acceptance of xxx level
    % -> largest qt_cutoff that gives us at least ar_cutoff > xxx
    [opt_cutoff, opt_qt] = cutoff_search(mat_Svec, smc_target_ar, smc_eff_zero);

    % --- (b) Correction
    mat_w = all(opt_cutoff<mat_Svec); %unnormalized weight
    mat_W = smc_nsim * (mat_w / sum(mat_w));

    % --- (c) Selection (we will do resampling all the time because some particles
    % are dead for sure)
    ESS = smc_nsim / (mean(mat_W.^2)); %should be roughly target_ar by construction

    % resampling (randomly select those that survived)
    ind_survived = find(all(opt_cutoff<mat_Svec));
    ind_selected = randsample(ind_survived, smc_nsim, true);
    mat_W = ones(1, smc_nsim); %reset as we resampled

    % --- (d) Mutation
    mat_B_old = mat_B; %copy
    mat_S_old = mat_S; %copy
    mat_R_old = mat_R; %copy
    mat_Q_old = mat_Q; %copy
    mat_X_old = mat_X; %copy

    parfor i=1:smc_nsim

        % current draws
        sind  = ind_selected(i);
        Bdraw = mat_B_old(:,:,sind);
        Sigmadraw = mat_S_old(:,:,sind);
        Rdraw = mat_R_old(:,sind);
        Qdraw = mat_Q_old(:,:,sind);
        Xdraw = mat_X_old(:,sind);

        % transition
        for j=1:smc_ntran
            [Bdraw,Sigmadraw,Rdraw,Qdraw,Xdraw,n_try_B, n_try_S, n_try_Q] = ess_sampler(Bdraw,Sigmadraw,Rdraw,Qdraw,Xdraw,function_restrictions_Svec, opt_cutoff, info, info_ess);
        end

        % check reponses that are restricted
        Svec = function_restrictions_Svec(Bdraw, Sigmadraw, Qdraw, info);

        % store result
        mat_B(:,:,i) = Bdraw;
        mat_S(:,:,i) = Sigmadraw;
        mat_R(:,i)   = Rdraw;
        mat_Q(:,:,i) = Qdraw;
        mat_X(:,i)   = Xdraw;
        mat_Svec(:,i) = Svec;
    end


    % --- stopping rule
    disp('===============');
    disp(['SMC Iteration = ', num2str(smc_iter)]);
    disp('opt_cutoff');
    disp(opt_cutoff);

    if smc_iter >= smc_max_iter
        disp('smc: max_iter reached');
        break
    end
    if abs(max(abs(opt_cutoff)) - smc_eff_zero) < (smc_eff_zero * 1e-2)
        disp('smc: successfully reached the desired posterior')
        break
    end

    % check if any of particles satisfy restrictions already
    ind_final = all(-smc_eff_zero < mat_Svec);
    if smc_early_stop && (sum(ind_final) > 2)
        disp('found a draw that satisfies restrictions');
        rst.mat_B = mat_B(:,:,ind_final);
        rst.mat_S = mat_S(:,:,ind_final);
        rst.mat_R = mat_R(:,ind_final);
        rst.mat_Q = mat_Q(:,:,ind_final);
        rst.mat_X = mat_X(:,ind_final);
        rst.mat_Svec = mat_Svec(:,ind_final);
        return
    end


    % Prep for the next iteration
    smc_iter = smc_iter + 1;
end

%% Collect all results
rst.mat_B = mat_B;
rst.mat_S = mat_S;
rst.mat_R = mat_R;
rst.mat_Q = mat_Q;
rst.mat_X = mat_X;
rst.mat_Svec = mat_Svec;
