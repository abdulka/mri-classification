%% [Kd,Yd,theta] = gprkerneldetrend(K,X,DATA,DT,Y,COVARFUNC) adjusts the kernel
%  matrix K for the covariets X, where the argument DATA=X*X' can be empty. The
%  adjustment is based on the index vector DT. Y is the raw input data
%  and COVARFUNC the optional covariance function, that is required to return the
%  covariance, given DATA, the gradient, and the Hession w.r.t. the
%  hyperparameters. If no covariance function is provided, the default is a
%  squared exponential with linear and constant component. Residual noise
%  parameter is always estimated as addiditonal hyperparameter.
%
%  Ahmed Abdulkadir, November 2014
%
%  $Id: gprkerneldetrend.m 298 2014-12-10 18:54:31Z ahmed.abdulkadir@cs.uni-freiburg.de $

% ----------------------------------------------------------------------- %
function [Kd,mu,s2,theta,optimout] = gprkerneldetrend(K, ... kernel matrix (feature space)
    X,                                       ... raw covariates
    data,                                    ... dot-product matrix of covariates
    dt,                                      ... indices of examples to use for parameter estimation
    Y,                                       ... raw data (input)
    covariancefunc)                          ... covarinace function

% set covariance function if not passed as argument
if ~exist('covariancefunc','var'), covariancefunc = @defaultcovariancefunc; end

% convert to cell is necessary
if ~iscell(covariancefunc), covariancefunc = {covariancefunc}; end

if size(K,1) == size(K,2) % if the first input is a square matrix

% --- BEIN OPTIMIZATION ------------------------------------------------------ %
% D: dimensionality of data
D = size(Y,2);

% theta0, lb, ub: initial guess of the hyperparameters and bounds
[theta,lb,ub]  = feval(covariancefunc{:},X(dt,:),data(dt,dt));
theta0         = [theta; 1];

% NOTE: fminunc and lbfgs give very similar results. L-BFGS-B was
% about 30% faster in most cases. Still, convergence was better with
% fminunc.
fun = @(theta) negloglikelihoodfunction(X,data,K,dt,D,theta,covariancefunc);

% **** INITIALIZE and RUN fminunc SOLVER *******
PlotFcns     = {@optimplotx @optimplotfunccount @optimplotfval ...
                @optimplotstepsize @optimplotfirstorderopt};
opts_fminunc = optimoptions(optimoptions('fminunc'), ...
    'Algorithm','trust-region', ...
    'DerivativeCheck','off',    ...
    'Diagnostics','off',        ...
    'DiffMaxChange',Inf,        ...
    'DiffMinChange',0,          ...
    'Display','iter',           ...
    'FinDiffRelStep',1e-8,      ...
    'FinDiffType','central',    ...
    'FunValCheck','on',         ...
    'GradObj','on',             ...
    'MaxFunEvals',5000,         ...
    'MaxIter',500,              ...
    'OutputFcn',[],             ...
    'PlotFcns',PlotFcns,        ...
    'TolFun',1e-6,              ...
    'TolX',1e-6,                ...
    ...'TypicalX',theta0,          ...
    ... trust-region Algorithm Only
    'Hessian','on',             ...
    ...'HessMult',[],           ...
    ...'HessPattern',[],        ...
    'MaxPCGIter',max(1,ceil(numel(theta0)/2)),             ...
    'PrecondBandWidth',5,     ...
    'TolPCG',0.1);

% this procudere will converge to a single local minimum
[theta,~,exitflag,output,grad,hessian] = fminunc(fun,gather(theta0),opts_fminunc);
assert(exitflag>0,'fminunc exitflag was: %i',exitflag);


% display final solution
%fprintf('theta=['); fprintf('%.2e ',theta); fprintf('\b]\n');

% C: covariance (without observation noise term)
C = feval(covariancefunc{:},X,data,theta,true(size(dt)));
% N: total number of examples
N = numel(dt);
% L: Cholesky decoposed kernel
[L,isnotpod] = chol(C(dt,dt) + theta(end)^2*eye(sum(dt)),'lower');
assert(~isnotpod,'Matrix is not positive definite.');
iC = (L^-1)'*L^-1;
% T: training indices in matrix form
T = full(sparse(1:numel(find(dt)),find(dt),ones(sum(dt),1)',sum(dt),N));
% R:  residual-forming matrix
R = (eye(N) - C(:,dt)*iC*T);
% Kd: detrended kernel
Kd = R*K*R';
% --- END OPTIMIZATION ------------------------------------------------------- %
else % otherwise assume that the first input represents the hyperparameters
    Kd    = [];
    optimout = [];
    theta = K;
    % C: covariance (without observation noise term)
    C = feval(covariancefunc{:},X,data,theta,true(size(dt)));
    % L: Cholesky decoposed kernel
    [L,isnotpod] = chol(C(dt,dt) + theta(end)^2*eye(sum(dt)),'lower');
    assert(~isnotpod,'Matrix is not positive definite.');
    iC = (L^-1)'*L^-1;
end


% compute posterior mean
if nargout > 1,
   mu = C(dt,:)'*(iC*Y(dt,:));
end
% compute posterior variance
if nargout > 2,
   % NOTE: the variance does not depend on the training or test targets,
   %       but is uniquely defined by the covariance
   s2 = bsxfun(@minus,diag(C),diag((C(:,dt)*iC*C(dt,:))));
end

% output optimization exit state
if (nargout > 4) && ~exist('optimout','var'),
    optimout = struct('exitflag',exitflag,'output',output,'grad',grad,...
        'hessian',hessian);
end
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% negative log marginal likelihood function, it's first and second 
% derivatives with respect to the hyperparameters theta, given a covariance
% function. 
function [nlik,g,H] = negloglikelihoodfunction(X,XX,tt,dt,D,theta,covariancefunc)

% inner product of training
tt = tt(dt,dt);

% 
if nargout > 2
    [C,dC,ddC] = feval(covariancefunc{:},X,XX,theta(1:end-1),dt);
    E          = eye(size(C),'like',C);
    dCdSigma   = 2*theta(end)*E;
    ddCddSigma = 2*E;
elseif nargout > 1
    [C,dC]     = feval(covariancefunc{:},X,XX,theta(1:end-1),dt);
    E          = eye(size(C),'like',C);
    dCdSigma   = 2*theta(end)*E;
else
    [C]        = feval(covariancefunc{:},X,XX,theta(1:end-1),dt);
    E          = eye(size(C),'like',C);
end

C = C + theta(end)^2*eye(size(C));

% N: number of examples
n = size(C,1);

% L: Cholesky decoposed covariance
% NOTE: Cholesky decomposition and subsequent inversion was about 10%
%       faster on the CPU compared with GPU (with a 360x360 matrix).
[L,isnotpd] = chol(C,'lower');
assert(~isnotpd, 'Matrix was not positive definite.');
if isnotpd
    % inverse of C
    e  = 1e-2;
    iC = inv(C+e*diag(n));
    warning('Regularizing with %e\ntheta_5=%e',e,theta(end));
    % log determinant
    logdet = 1/2*log(det(C));
else
    % inverse of lower triangula part
    iL = inv(L); % this line takes almost one third of the function time
    % iC: inverse of covariance matrix
    iC = iL'*iL;
    % log determinant
    logdet = sum(log(diag(L)));
end
if ~isreal(logdet)
    assert(isreal(logdet),'Assertion failed! Determinant is not rela %f+%fi',     ...
        real(logdet),imag(logdet));
    logdet = real(logdet);
end

% lik: likelihood
% complexity + constant + data
nlik = double(gather(D*logdet + n*D/2*log(2*pi) + 1/2*sum(sum(iC.*tt))));

if nargout > 1
    % g: gradient
    if isa(dC{1},'gpuArray')
        g = gpuArray.zeros(numel(dC+1),1,classUnderlying(dC{1}));
        for i=1:numel(dC)
            g(i) = +D/2*sum(sum(iC.*dC{i})) ...
                -1/2*sum(sum((iC*dC{i}*iC).*tt));
        end
        g(end) = +D/2*sum(sum(iC.*dCdSigma)) ...
                -1/2*sum(sum((iC*dCdSigma*iC).*tt));
        % push to CPU memory
        g = gather(g);
    else % we assume that it's on the CPU
        % g: gradient
        g = cellfun(@(dC) +D/2*sum(sum(iC.*dC)) ...
            -1/2*sum(sum((iC*dC*iC).*tt)), [dC{:} {dCdSigma}], ...
            'UniformOutput',true)';
    end
end

if nargout > 2
    % H: Hassian matrix
    if isa(dC{1},'gpuArray')
        H = gpuArray.zeros(numel(theta),classUnderlying(tt));
        for i=1:numel(dC)
            for j=1:i
                H(i,j) = +D/2*sum(sum(iC.*ddC{i,j} - (iC*dC{j}*iC).*dC{i}))      ...
                    -1/2*sum(sum(( iC*(-dC{j}*iC*dC{i} + ddC{i,j} ...
                    -dC{i}*(iC*dC{j}))*iC).*tt));
                H(j,i) = H(i,j);
            end
        end
        H = gather(H);
    else
        % H: Hessian matrix, initialized at zeros
        H      = zeros(numel(theta),'like',tt);
        
        % second derivatives with respect to the hyperparameters
        % NOTE: it should be posible to make this more efficient by
        % factorizing correctly
        for i=1:numel(dC)
            for j=1:i
                H(i,j) = +D/2*sum(sum(iC.*ddC{i,j} - (iC*dC{j}*iC).*dC{i})) ...
                    -1/2*sum(sum((iC*(-dC{j}*iC*dC{i} + ddC{i,j} ...
                    -dC{i}*(iC*dC{j}))*iC).*tt));
                H(j,i) = H(i,j);
            end
        end
        H(end,end) = +D/2*sum(sum(iC.*ddCddSigma - (iC*dCdSigma*iC).*dCdSigma)) ...
                    -1/2*sum(sum((iC*(-dCdSigma*iC*dCdSigma + ddCddSigma ...
                    -dCdSigma*(iC*dCdSigma))*iC).*tt));
    end
end
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
%% Squared exponential covariance function with linear component
%  (default covariance function)
%
%  Inputs:
%  data     linear kernel matrix
%  theta    vector of M hyperparameters
%  idx      train indices
%
%
%  Initialization and default bounds on parameters theta:
%  [theta0 lb ub] = squaredexponentialcovariance(data)
%  theta0   stable initial paramaters
%  lb       lower bound for theta
%  ub       upper bound
%
%  Evaluation of covariance function:
%  [K,C,dC,ddC] = squaredexponentialcovariance(data,theta,idx)
%  K        evaluated kernel
%  C        evaluated kernel with measurement noise (theta(4))
%  dC       cell array of derivatives of the covariance with respect to the
%           hyperparameters theta
%  ddC      cell array of second derivatives of the covariance with respect
%           to the hyperparameters theta

%  idx:      training examples
function [C,dC,ddC] = defaultcovariancefunc(~,data,theta,idx)

% n: number of parameters
n  = 4;

% return
if nargin == 0
    C = n;
    return;
end

% assert that X is a numeric array, i.e. linear kernel
assert(isnumeric(data),'Error: ''data'' must be a numeric array.');

% theta0: initialize parameters and bounds
if nargin < 3
    C  = +ones(n,1, 'like',data); % initial parameters
    % estimate reasonable RBF width parameter where the mean of the off-diagonal
    % elements of the kernel are close to 0.5
    m    = size(data,1);
    D    = diag(data)*ones(1,m) -2*data  + ones(m,1)*diag(data)';
    kernelwidths = exp(linspace(-16,16,200));
    meanoffdiag = arrayfun(@(theta) (sum(sum(exp(-1./(theta^2).*D)))-m)/(m*m-m),  ...
        kernelwidths);
    [~,imin] = min(abs(meanoffdiag-0.25));
    C(2) = kernelwidths(imin);
    C(3) = sqrt(1/mean(abs(diag(data)))); 
    C(4) = 0;
    C(1) = 1;
    
    dC  = -Inf(n,1,'like', data); % lower bound
    ddC = +Inf(n,1,'like', data); % upper bound
    return;
end

% idx: training indices
if exist('idx','var')
    data = data(idx,idx);
end

% m: number of examples
m = size(data,1);

% D: matrix of pair-wise squared Euclidean distances
D = diag(data)*ones(1,m) -2*data  + ones(m,1)*diag(data)';

% C: covariance
B = exp(-1./(theta(2)^2).*D);
C = theta(1)^2.*B + theta(3)^2.*data + theta(4)^2;

% C: iid noise-augmented covariance
if nargout > 1
    Z = zeros(m,'like',data);
    O = ones( m,'like',data);
end

% dC: derivative of C w.r.t. theta
if nargout > 1
    dC = { 2*theta(1)*B                    ,  ... scale of RBF portion
           2*theta(1)^2/(theta(2)^3).*D.*B ,  ... width of RBF
           2*theta(3)*data                 ,  ... scale of linear portion
           2*theta(4)*O                    }; ... offset
     % noise
end

% ddC: second derivative of C w.r.t theta
if nargout > 2
    ddC = repmat({ Z },n,n);
    ddC{1,1} = 2*B;
    ddC{1,2} = 4*theta(1)/(theta(2)^3).*D.*B;
    ddC{2,1} = ddC{1,2};
    ddC{2,2} = 2*theta(1)^2.*B.*D.*(2*D-3*theta(1)^2)/(theta(2)^6);
    ddC{3,3} = 2*data;
    ddC{4,4} = 2*O;
end
% ----------------------------------------------------------------------- %
