function density = inverse_wishart_pdf(X, Psi, nu)
% INVERSE_WISHART_PDF  Evaluates the probability density function of an Inverse Wishart distribution
%
% INPUTS:
%   X   - p×p symmetric positive definite matrix at which to evaluate the density
%   Psi - p×p symmetric positive definite scale matrix
%   nu  - degrees of freedom (nu > p-1)
%
% OUTPUT:
%   density - value of the Inverse Wishart density at X
%
% The Inverse Wishart density is given by:
%   f(X) = |Psi|^(nu/2) |X|^(-(nu+p+1)/2) exp(-0.5*trace(Psi*inv(X))) / (2^(nu*p/2) * pi^(p(p-1)/4) * prod_{i=1}^p Gamma((nu+1-i)/2))
%
% where p is the dimension of X.

    % Check if X is square
    [p, q] = size(X);
    if p ~= q
        error('Input matrix X must be square');
    end
    
    % Check if X is symmetric
    if norm(X - X', 'fro') > 1e-10
        error('Input matrix X must be symmetric');
    end
    
    % Check if nu is valid
    if nu <= p - 1
        error('Degrees of freedom nu must be greater than p-1');
    end
    
    % Compute the normalizing constant
    log_norm_const = (nu*p/2)*log(2) + (p*(p-1)/4)*log(pi);
    for i = 1:p
        log_norm_const = log_norm_const + log(gamma((nu+1-i)/2));
    end
    
    % Compute the log-density
    log_det_Psi = log(det(Psi));
    log_det_X = log(det(X));
    trace_term = trace(Psi / X);
    
    log_density = (nu/2)*log_det_Psi - ((nu+p+1)/2)*log_det_X - 0.5*trace_term - log_norm_const;
    
    % Convert log-density to density
    density = exp(log_density);
end