function [yf,exogf,yf_unc,e_tplus1_tplush,KL2] = DrawConditionalForecast(y,A0,B_draw,exog,yCondition,eCondition,Omega_f,Omega_g)
%% **************** FUNCTION: DRAWCONDITIONALFORECAST ********************%
% This function draws a conditional forcast from the truncated distribution
% according the variables conditioned upon. The function encompassed three
% types of conditioning: Conditioning on a path of one or more observables;
% Conditioning on a path of one or more shocks; or Conditioning on a path
% of both shocks and obserables.
%
% Inputs:   o y:            The data
%           o A0:           The orthogonalizing matrix from the sVAR
%           o B_draw:       A draw from the coefficients coming from the
%                           reduced-form BVAR
%           o exog:         The constant in the BVAR       
%           o yCondition:   A pxn matrix of the path for observables (p =
%                           no. of restriction; n = no. of observables.
%                           Columns for observables that one does not want
%                           to condition on require a column of NaNs.
%           o eCondition:   A pxn matrix of the path for shocks (p =
%                           no. of restriction; n = no. of observables
%                           Columns for shocks that one does not want
%                           to condition on require a column of NaNs.
%           o Omega_f:      VCov for the conditioned variable
%           o Omega_g:      VCov for the conditioned shock (note if one
%                           wants to condition-on-structural-scenario,
%                           this must be of size n-p x n-p)
%
% Outputs:  o yf:           Conditional forecast of the observables
%           o exogf:        constant
%           o yf_unc:       unconditional forecast
%           o prob_cond:    likelihood of conditional forecast from
%                           unconditional distribution
%           o prob_unc:     likelihood of unconditional forecast from
%                           unconditional distribution
% ------------------------------------------------------------------------%          
%% Find type of conditional forecast

if ~isempty(yCondition) && isempty(eCondition)
    type = 'condition-on-observables';
    h = size(yCondition,1);
elseif isempty(yCondition) && ~isempty(eCondition)
    type = 'condition-on-shock';   
    h = size(eCondition,1);
elseif ~isempty(yCondition) && ~isempty(eCondition)
    type = 'condition-on-structural-scenario';
    h = size(yCondition,1);
else
    error('No conditions inputted. Please insert condition paths into yCondition and/or eCondition')
end

%% Calculate objects from unconditional forecasts using compact notation
%(Section 1.2 - Unconditional Forecasting)

[T,n] = size(y);
[b_tplus1_tplush,M] = compactForecastNotation(A0,B_draw,exog,T,h,y);

%% Calculate objects depending on conditions

switch type
    case 'condition-on-observables' % (Section 1.3)

        y_tplus1_tplush = vec(yCondition');
        C = diag(isfinite(y_tplus1_tplush));
        C(isnan(y_tplus1_tplush),:) = [];
        
        D = C*M';
        if isempty(Omega_f) 
            Omega_f = D*D';
        end
        
        f_tplus1_tplush = y_tplus1_tplush(isfinite(y_tplus1_tplush));
        Omega_f_final = Omega_f;
        
    case 'condition-on-shock' % (Section 1.4)
        e_tplus1_tplush = vec(eCondition');
        Psi = diag(isfinite(e_tplus1_tplush));
        Psi(isnan(e_tplus1_tplush),:) = [];
        C = Psi/M';
        D = C*M';
        
        e_tplus1_tplush(isnan(e_tplus1_tplush))=0;
        f_tplus1_tplush = C*b_tplus1_tplush + Psi*e_tplus1_tplush;
        Omega_f_final = Omega_g;
        
    case 'condition-on-structural-scenario'  % (Section 1.5)
        
        y_tplus1_tplush = vec(yCondition');
        f_tplus1_tplushY = y_tplus1_tplush(isfinite(y_tplus1_tplush));
        C_upper_bar = diag(isfinite(y_tplus1_tplush));
        C_upper_bar(isnan(y_tplus1_tplush),:) = [];
        
        D_upper_bar = C_upper_bar*M';
        if isempty(Omega_f)
            Omega_f = D_upper_bar*D_upper_bar';
        end
        
        
        e_tplus1_tplush = vec(eCondition');
        Psi = diag(isfinite(e_tplus1_tplush));
        Psi(isnan(e_tplus1_tplush),:) = [];
        
        C_lower_bar = Psi/(M');
        C = [C_upper_bar ; C_lower_bar];
        D = C*M';
        
        f_tplus1_tplush = [f_tplus1_tplushY ; C_lower_bar*b_tplus1_tplush];
        Omega_f_final = blkdiag(Omega_f,Omega_g);
        
end

%% Calculate Mean and Variance of conditional and unconditional forecasts



Dstar = pinv(D);  % Dstar is the generalised inverse of D
Dhat = null(D)';  % The rows of Dhat form an orthonormal basis for the null space of D

mu_y_unc = b_tplus1_tplush;
Omega_y_unc = (M'*M);

mu_y = M'*Dstar*f_tplus1_tplush+M'*(Dhat'*Dhat)/(M')*b_tplus1_tplush; % Equation (10)
Omega_y = M'*(Dstar*Omega_f_final*Dstar'+(Dhat'*Dhat))*M; % Equation (10)
Omega_y = (Omega_y+Omega_y')*.5; 

mu_e = Dstar*f_tplus1_tplush - Dstar*C*b_tplus1_tplush;
Omega_e = Dstar*Omega_f_final*Dstar'+Dhat'*Dhat;
Omega_e = (Omega_e+Omega_e')*.5; 

%% Draw objects from conditional Posterior Distribution 
y_tplus1_tplush = mvnrnd(mu_y,Omega_y)';
y_tplus1_tplush_unc = mvnrnd(mu_y_unc,Omega_y_unc)';

e_tplus1_tplush = ((y_tplus1_tplush' - b_tplus1_tplush')/(M))';
e_tplus1_tplush = reshape(e_tplus1_tplush,[n,h])';
y_tplus1_tplush = reshape(y_tplus1_tplush,[n,h])';

y_tplus1_tplush_unc = reshape(y_tplus1_tplush_unc,[n,h])';

yf = [y; y_tplus1_tplush];
yf_unc = [y; y_tplus1_tplush_unc];

exogf = [exog; ones(size(y_tplus1_tplush,1),1)];

%% Kullback Leibler Divergence

% SigmaEpsInv = inv(Omega_e);

KL2 = 0.5*(trace(Omega_e)+mu_e'*mu_e-n*h+log(1/det(Omega_e)));

% KL_calibrated = (1 + (1 - exp(-2*mean(KL_save)/(n*hmax)))^(0.5))/2;

end

