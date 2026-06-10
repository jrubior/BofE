function [y, betahat, sigmahat] = negative_loglike_obj_mp_with_alpha(x,info,Xnd,Ynd,MIN,MAX,SS,covid_periods,Tcovid)

  T = size(Ynd,1);
  d= info.nvar +2; 
  idx = 1;

  % llambda (siempre primero)
  llambda0 = MIN.llambda + (MAX.llambda - MIN.llambda) ./ (1 + exp(-x(idx)));
  idx = idx + 1;

  % ahora vienen o bien ppsi (GLP, length = nvar) o eta (LP, length = covid_periods+1)
    if isempty(Tcovid)
        % GLP: ppsi (one per variable)
        ppsi_idx = idx : (idx + info.nvar - 1);
        ppsi0 = MIN.ppsi + (MAX.ppsi - MIN.ppsi) ./ (1 + exp(-x(ppsi_idx)));
        idx = idx + info.nvar;
    else
        % LP: eta vector (covid_periods + 1)
        ncp = covid_periods + 1;
        eta_idx = idx : (idx + ncp - 1);
        eta = MIN.eta' + (MAX.eta' - MIN.eta') ./ (1 + exp(-x(eta_idx)));
        idx = idx + ncp;
    end

    alpha0 = MIN.alpha + (MAX.alpha - MIN.alpha) ./ (1 + exp(-x(idx)));
    idx = idx + 1;
    
    if info.use_soc == 1
        miu0 = MIN.miu + (MAX.miu - MIN.miu) ./ (1 + exp(-x(idx)));
        idx = idx + 1;
    end

    if info.use_dio == 1
        theta0 = MIN.theta + (MAX.theta - MIN.theta) ./ (1 + exp(-x(idx)));
        idx = idx + 1;
    end
 

  if ~isempty(Tcovid)
      ppsi0  = SS; 
      
      invweights = ones(T,1);

      covid_periods_eff = min(covid_periods, T - Tcovid + 1);
      if covid_periods_eff < covid_periods
        warning('covid_periods truncated to fit sample length.');
      end
        
     covid_idx = Tcovid : (Tcovid + covid_periods_eff - 1);
     invweights(covid_idx) = eta(1:covid_periods_eff);
        
        rho = eta(end);                        % rho siempre el último
        post_start = Tcovid + covid_periods_eff;
        
        if post_start <= T
            eta_last = eta(covid_periods_eff);
            h = (1:(T - post_start + 1))';     % empieza en 1 (NO repite eta_last)
            invweights(post_start:T) = 1 + (eta_last - 1) .* (rho .^ h);
        end
        
        Ynd = diag(1./invweights) * Ynd;
        Xnd = diag(1./invweights) * Xnd;


  end
  

  % Prior matrices and vectors
  nnuBar              = info.nvar + 2;
  PpsiBar             = diag(ppsi0);
  oomegaBar = zeros(info.m, 1);
  oomegaBar(1:info.nex,1) = info.Vc;
  
  for ell=1:info.nlag
    oomegaBar(info.nex + 1 + (ell - 1)*info.nvar:info.nex + ell*info.nvar, 1) =  (llambda0^2)*(nnuBar-info.nvar-1)./((ell^alpha0)*ppsi0); % added alpha
  end
      
  OomegaBar=diag(oomegaBar);
  OomegaBarInverse=diag(1./oomegaBar);
  mmuBar = zeros(info.m,info.nvar);
  diagmmuBar=ones(info.nvar,1);
  diagmmuBar(info.pos)=0;   % Set to zero the prior mean on the first own lag for variables selected in the vector pos
  mmuBar(info.nex+1:info.nvar+info.nex,:)=diag(diagmmuBar); 
  
  % compute the marginal likelihood
  Y = Ynd;
  X = Xnd;
  
  Td= 0;
  ydsoc= [];
  xdsoc = [];
  yddio = [];
  xddio = [];

  if info.use_soc == 1 
      scl = 1/miu0;                                   % convención: (1/miu)
      ydsoc = scl * diag(info.y0bar);               % n x n diagonal
      if isfield(info,'pos') && ~isempty(info.pos)
          ydsoc(info.pos, info.pos) = 0;            % respetar pos si aplica
      end
      xdsoc = [zeros(info.nvar,1), scl*repmat(diag(info.y0bar),1,info.nlag)]; % n x (nex + nvar*nlag)
      Y = [Y; ydsoc];
      X = [X; xdsoc];
      Td = Td + size(Y,2);
      
  end
  
   if info.use_dio == 1
    
    scl = 1/theta0;   % cuidado: theta0 debe estar definido antes
    yddio = (scl * info.y0bar(:))';   % forzamos fila
    if isfield(info,'pos') && ~isempty(info.pos)
        yddio(info.pos) = 0;            % respetar pos si aplica
    end
    % xddio: 1 x (1 + nvar*nlags)  (constante + repeticiones de y0)
    xddio = [scl, scl * repmat(info.y0bar(:)', 1, info.nlag)];

    % anexar exactamente como hace sur
    Y = [Y; yddio];
    X = [X; xddio];

    % comportamiento idéntico al if sur...: asignar Td = 1
    Td = Td + 1;
  end



  T= T+Td; 

  betahat = (X'*X  + OomegaBarInverse)\(X'*Y + OomegaBarInverse*mmuBar);
  epshat = Y - X*betahat;  
  S = epshat'*epshat + PpsiBar + (betahat - mmuBar)'*OomegaBarInverse*(betahat - mmuBar);
  sigmahat = S / (T + d + info.nvar + 1);

  
  % numerical stability tricks
  DPsi = chol(PpsiBar\eye(info.nvar),'lower'); %diag(1./sqrt(psi)), en logformin
  D = chol(OomegaBar,'lower'); %diag(sqrt(omega))
  
  aaa = D'*(X'*X)*D;
  bbb = DPsi'*(epshat'*epshat  + (betahat - mmuBar)'*OomegaBarInverse*(betahat - mmuBar))*DPsi;
  
  eigaaa = real(eig(aaa));
  eigaaa(eigaaa < 1e-12) = 0;
  eigaaa = eigaaa + 1;

  eigbbb = real(eig(bbb)); 
  eigbbb(eigbbb<1e-12) = 0; 
  eigbbb = eigbbb + 1;


  tmp1_stable = - T*0.5*LogAbsDet(PpsiBar) - (T+nnuBar)*0.5*sum(log(eigbbb));
  tmp2_stable = - info.nvar*0.5*sum(log(eigaaa));
        
  tmpg = 0;
  for i = 1:info.nvar
    tmpg =  tmpg  + gammaln((T + nnuBar + 1 - i)/2) - gammaln((nnuBar + 1 - i)/2);
  end
        
  yextended_stable = - info.nvar*T*0.5*log(pi) + tmpg + tmp1_stable + tmp2_stable; %+ - T*sum(log(ppsi0*(nnuBar-info.nvar-1)))/2;
            
  y = (yextended_stable);

  if info.use_soc == 1 || info.use_dio ==1
       yd=[ydsoc;yddio];
       xd=[xdsoc;xddio];

      betahatd= mmuBar;
      epshatd = yd- xd * betahatd; 
      aaa = D'*(xd'*xd)*D;
      bbb = DPsi'*(epshatd'*epshatd  + (betahatd - mmuBar)'*OomegaBarInverse*(betahatd - mmuBar))*DPsi;
      eigaaa=real(eig(aaa)); eigaaa(eigaaa<1e-12)=0; eigaaa=eigaaa+1;
      eigbbb=real(eig(bbb)); eigbbb(eigbbb<1e-12)=0; eigbbb=eigbbb+1;
      % normalizing constant
      norm = - info.nvar*Td*log(pi)/2 + sum(gammaln((Td+d-[0:info.nvar-1])/2)-gammaln((d-[0:info.nvar-1])/2)) +...
           - Td*sum(log(ppsi0))/2 - info.nvar*sum(log(eigaaa))/2 - (T+d)*sum(log(eigbbb))/2;
      y = y-norm; 
  end 

  if ~isempty(Tcovid)
        y = y - info.nvar*sum(log(invweights));
  end 

  % Priors for hyperparameters 
  y = y + logGammapdf(llambda0, info.priorcoef.llambda.k, info.priorcoef.llambda.theta);
  
  if ~isempty(Tcovid)
    y = y - 2*sum(log(eta(1:covid_periods_eff))) + logBetapdf(eta(covid_periods_eff+1), info.priorcoef.eta4.alpha, info.priorcoef.eta4.beta);
  end
  if isempty(Tcovid)
    y = y + sum(logIG2pdf(ppsi0/(nnuBar-info.nvar-1), info.priorcoef.alpha.PSI, info.priorcoef.beta.PSI));
  end
  
  if info.use_soc == 1
      y= y +logGammapdf(miu0,info.priorcoef.miu.k,info.priorcoef.miu.theta);
  end
   if info.use_dio == 1
      y= y +logGammapdf(theta0,info.priorcoef.theta.k,info.priorcoef.theta.theta);
  end

  y = y*-1;


function r=logGammapdf(x,k,theta);
r=(k-1)*log(x)-x/theta-k*log(theta)-gammaln(k);

function r=logIG2pdf(x,alpha,beta);
r=alpha*log(beta)-(alpha+1)*log(x)-beta./x-gammaln(alpha);

function r=logBetapdf(x,al,bet);
r=(al-1)*log(x)+(bet-1)*log(1-x)-betaln(al,bet);