%% Kd = kernelLinearDetrend(K,X,IDX,LMBD) linearly detrends kernel K using the
%  regressors in X. The logical array IDX indicates training set used to
%  estimate the parameters for the detrend. The parameter LMBD is the ridge
%  parameter that biases the solution towards zero in order to reduce
%  overfitting.
%
%  Usage Example:
%    X = [ones(100,1) [linspace(0,20,40)'; 20*rand(20,1) ;linspace(0,20,40)']];
%    Y = [X*[1 2]' [ones(50,8); -ones(50,8)] 0.6*randn(100,10)];
%    K = Y*Y';
%    idx = true(100,1);
%    idx(41:60) = false;
%    lmbd = 0.01; % ridge parameter
%    Kd = kernelLinearDetrend(K,X,idx,lmbd);
%    figure(1);
%    subplot(1,2,1); imagesc(K);axis image;
%    subplot(1,2,2); imagesc(Kd); axis image
%
%  Ahmed Abdulkadir, February 2013
%
% http://www.m-hikari.com/ams/ams-2010/ams-9-12-2010/dorugadeAMS9-12-2010.pdf

function [Kd,R] = kernelLinearDetrend(K,X,idx,lmbd)

% lmbd: ridge parameter
if ~exist('lmbd','var') || isempty(lmbd); lmbd = 0; end

% X: detrending matrix
% !! the 'gather' function is required because 'pinv.m' is not compatible
% with gpuArrays (as of R2013a)
X = gather(X);

% N: number of exapmles
N = size(K, 1);

% check if idx is defined and take all examples as training otherwise
if ~exist('idx','var') || isempty(idx),
    idx = true(N,1);
else
    % idx: assure logical array
    idx = logical(idx);
end

% assert compatible dimensions
assert(size(K, 2) == N, 'Kernel matrix K must be square.');
assert(size(X, 1) == N, 'size(X,1) must be equal to size(K,1).');
assert(numel(idx) == N, 'numel(idx) must be equal to size(K,1).');

% training kernel including all examples
if nargout > 1
    [Kd,R] = doKernelDetrending(K,X,idx,lmbd);
else
    Kd = doKernelDetrending(K,X,idx,lmbd);
end

end

function [Kd,R] = doKernelDetrending(K,X,idx,lmbd)

% N: number of elements
N = numel(idx);
% n: number of example used to estimate parameters
n = sum(idx);
% p: number of columns in X
p = size(X,2);

% T: training indices in matrix form in order to use a restricted set of
%    examples for estimation of the regression parameters.
T = full(sparse(1:n,find(idx),ones(n,1)',n,N));

% R:  residual-forming matrix
if lmbd>0,
    pinvX = (X(idx,:)'*X(idx,:) + lmbd*eye(p))^-1*X(idx,:)';
else
    pinvX = pinv(gather(X(idx,:)));
end
R = (eye(N) - X*pinvX*T);

% compute detrended kernel
Kd = R*K*R';

end

