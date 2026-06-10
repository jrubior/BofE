function Q = draw_Q_conditional_on_sign(Sigmadraw, info)

% This function draws Q from the uniform conditional on sign restrictions and reduced form parameter
% for simplicity we limit the sign restrictions to be a function of Sigmadraw
% but, as you can see below, extension is fairly standard

% number of shocks identified
n_shock_identified = numel(info.Ss);

% covariance matrix
chol_lower_Sigmadraw = chol(Sigmadraw, 'lower');

% --- Q drawing
% idea: sequentially draw the column conditional on previous columns 
% draw new column from Haar in the orthogonal complement of the previous columns

Q = zeros(info.nvar, info.nvar);
% tic
for j = 1:info.nvar
    
    valid = 0;
    while(~valid)
        
        v = randn(info.nvar, 1);               % step 1: fresh Gaussian draw
        
        % sequential drawing
        if j > 1                               % step 2: remove previous columns
            v = v - Q(:,1:j-1) * (Q(:,1:j-1).' * v);
        end
        nrm = norm(v);                         % safety check
        if nrm < 1e-11
            error('Numerical rank deficiency—draw again');
        end
        v = v / nrm;                           % step 3: normalise
        if rand > 0.5,  v = -v;  end           % random sign flip (optional)
        Q_j = v;
        
        % step 4: check sign restrictions j-th shock
        if j <= n_shock_identified
            L0_j = chol_lower_Sigmadraw*Q_j;
            valid = all(info.Ss{j,1}*L0_j>0); % check sign restriction
        else
            valid = 1; % we only identify n_shock_identified
        end
    end
    Q(:,j) = Q_j;
end
% toc

%% Check if this is true
% S_prop = function_restrictions_i(Bdraw, Sigmadraw, Q,fo_inv,fo_str2irfs,info);


