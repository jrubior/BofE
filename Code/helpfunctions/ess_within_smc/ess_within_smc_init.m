function rst = ess_within_smc_init(data_num, function_restrictions_Svec, model, init_params)
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
reoptimizeMH = 1;

% Para reproducibilidad:
seed = 123; 

rng(seed,'twister'); 

% Cliente
mainStream = RandStream('Threefry','Seed',seed);
RandStream.setGlobalStream(mainStream);

% Workers
psc = parallel.pool.Constant(@() RandStream('Threefry','Seed',seed)); % Crea en cada worker su propio RandStream usando una fábrica

parfor_call_id = 0;
smc_nsim = 2000;
smc_eff_zero = 0; 
smc_target_ar =0.2;
smc_ntran=10;
smc_max_iter=1000;
smc_early_stop = 1;
%% Actual loop

% --- iid draw from the posterior (sitting here for testing and to get dimension)
cholPpsiTilde = chol(model.PpsiTilde,'lower');
inv_cholPpsiTilde_prime = inv(cholPpsiTilde');
cholOomegaTilde = chol(model.OomegaTilde,'lower'); % this matrix is used to draw B|Sigma below

[Bdraw, Sigmadraw, Rdraw, Qdraw, Xdraw] = iid_sampler(model.mmuTilde, cholOomegaTilde, inv_cholPpsiTilde_prime, model.nnuTilde, model);
Svec = function_restrictions_Svec(Bdraw, Sigmadraw, Qdraw, model);

% --- Setting up for ESS
% definition to facilitate the draws from B|Sigma
hh              = @(x)chol(x);
cholOomegaTilde = hh(model.OomegaTilde)'; % this matrix is used to draw B|Sigma below

chol_lower_OomegaTilde = chol(model.OomegaTilde, 'lower');
chol_lower_OomegaTilde_inv = eye(size(model.OomegaTilde))/chol_lower_OomegaTilde;

Rmean = nan(model.nvar,model.nnuTilde);
for i=1:model.nnuTilde
    Rmean(:,i) = zeros(model.nvar,1);
end
Rvariance = zeros(model.nvar*model.nnuTilde,model.nvar*model.nnuTilde);
for i=1:model.nnuTilde
    Rvariance((i-1)*model.nvar+1:i*model.nvar,(i-1)*model.nvar+1:i*model.nvar) = model.PpsiTilde\eye(model.nvar);
end
%chol_lower_Rvariance = chol(Rvariance,'lower');
chol_lower_Rvariance_common = chol(model.PpsiTilde\eye(model.nvar), 'lower');

% collect information for ess
info_ess = [];
info_ess.Rmean = Rmean;
info_ess.chol_lower_Rvariance_common = chol_lower_Rvariance_common;
info_ess.nnuTilde = model.nnuTilde;

info_ess.mmuTilde = model.mmuTilde;
info_ess.chol_lower_OomegaTilde = chol_lower_OomegaTilde;
info_ess.chol_lower_OomegaTilde_inv = chol_lower_OomegaTilde_inv;
info_ess.cholOomegaTilde = cholOomegaTilde;

% --- Initialization
mat_B = nan([model.pred,model.nvar, smc_nsim]);
mat_S = nan([model.nvar,model.nvar, smc_nsim]);
mat_R = nan([size(Rdraw) smc_nsim]);
mat_Q = nan([model.nvar,model.nvar, smc_nsim]);
mat_X = nan([size(Xdraw), smc_nsim]);

mat_Svec = nan([numel(Svec), smc_nsim]);
mat_W = ones(1, smc_nsim); % initial weihgting is 1 because they came from prior

for i=1:1:smc_nsim

    % iid draws
    [Bdraw, Sigmadraw, Rdraw, Qdraw, Xdraw] = iid_sampler(model.mmuTilde, cholOomegaTilde, inv_cholPpsiTilde_prime, model.nnuTilde, model);
    Svec = function_restrictions_Svec(Bdraw, Sigmadraw, Qdraw, model);

    % store
    mat_B(:,:,i) = Bdraw;
    mat_S(:,:,i) = Sigmadraw;
    mat_R(:,i)   = Rdraw;
    mat_Q(:,:,i) = Qdraw;
    mat_X(:,i)   = Xdraw;
    mat_Svec(:,i) = Svec;

end

if ~isempty(init_params.Q_init)

    Svec = function_restrictions_Svec(init_params.B_init, init_params.S_init, init_params.Q_init, model);

    % store
    mat_B(:,:,1) = init_params.B_init;
    mat_S(:,:,1) = init_params.S_init;
    mat_R(:,1)   = init_params.R_init;
    mat_Q(:,:,1) = init_params.Q_init;
    mat_X(:,1)   = init_params.X_init;
    mat_Svec(:,1) = Svec;

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

    parfor_call_id = parfor_call_id + 1;

    parfor i=1:smc_nsim
        % Para reproducibilidad con 'parfor' (según chatGPT)
        stream = psc.Value; % cada worker obtiene el stream compartido
        stream.Substream = (parfor_call_id - 1) * smc_nsim + i;
        RandStream.setGlobalStream(stream);

        % current draws
        sind  = ind_selected(i);
        Bdraw = mat_B_old(:,:,sind);
        Sigmadraw = mat_S_old(:,:,sind);
        Rdraw = mat_R_old(:,sind);
        Qdraw = mat_Q_old(:,:,sind);
        Xdraw = mat_X_old(:,sind);

        % transition
        for j=1:smc_ntran
            [Bdraw,Sigmadraw,Rdraw,Qdraw,Xdraw,n_try_B, n_try_S, n_try_Q] = ess_sampler(Bdraw,Sigmadraw,Rdraw,Qdraw,Xdraw,function_restrictions_Svec, opt_cutoff, model, info_ess);
        end

        % check reponses that are restricted
        Svec = function_restrictions_Svec(Bdraw, Sigmadraw, Qdraw,model);

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
    % if abs(max(abs(opt_cutoff)) - smc_eff_zero) < (smc_eff_zero * 1e-2)
    %     disp('smc: successfully reached the desired posterior')
    %     break
    % end

    % check if any of particles satisfy restrictions already
    ind_final = all(smc_eff_zero < mat_Svec);
    if smc_early_stop && (sum(ind_final) > 0) % Porque 2 y no 0?
        disp('found a draw that satisfies restrictions');
        rst.mat_B = mat_B(:,:,ind_final);
        rst.mat_S = mat_S(:,:,ind_final);
        rst.mat_R = mat_R(:,ind_final);
        rst.mat_Q = mat_Q(:,:,ind_final);
        rst.mat_X = mat_X(:,ind_final);
        rst.mat_Svec = mat_Svec(:,ind_final);
        rst.ind_final=ind_final;
        return
    end


    % Prep for the next iteration
    smc_iter = smc_iter + 1
end

%% Collect all results
rst.mat_B = mat_B;
rst.mat_S = mat_S;
rst.mat_R = mat_R;
rst.mat_Q = mat_Q;
rst.mat_X = mat_X;
rst.mat_Svec = mat_Svec;
rst.ind_final=ind_final;
