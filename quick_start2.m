% Use partial discrete Walsh-Hadamart transform matrix
% (see function Utilities/pdwht_operator.m)
%
% Use DCT as sparsifying basis 
% (requires signal processing toolbox)

clear; clc;
if ~exist('Utilities','dir'); Run_Me_1st; end
if ~exist('fWHtrans','file');
    error('This demo requires the discrete Walsh-Hadamart transform. Please compile fWHtrans.cpp in the folder Utilities by calling ''mex -O fWHtrans.cpp'''); 
end
rng('default');

% problem sizes 
n = 1024*8;     % must be a power of 2
m = n/8;
fprintf('[n,m] = [%i,%i]\n',n,m);

% generate xs (non-sparse)
xs = 100*cumsum(randn(n,1)); 

% A = partial DWHT matrix
p = randperm(n);
picks = sort(p(1:m),'ascend'); picks(1) = 1;
perm = randperm(n);
A = pdwht_operator(picks, perm);

% b = A*xs + noise
sigma = 0.2;    % noise level
noise = sigma*randn(m,1);
bs =  A.times(xs);
b = bs + noise;

% set options
tol = 5e-4;
rho = 5e-4;
opts.tol = tol;
opts.rho = rho;
opts.basis.times = @(x)  dct(x);    % sparsfying basis
opts.basis.trans = @(x) idct(x);    % inverse of basis

% call YALL1
tic; [x,Out] = yall1(A, b, opts); 
fprintf('\n');
toc
relerr = norm(x-xs)/norm(xs);
fprintf('iter = %4i, error = %e\n',Out.iter,relerr)
plot(1:n,xs,'b-',1:n,x,'r:'); 
legend('Original','Recovered','location','best')


save('data/quick_start2.mat', 'n', 'm', 'xs', 'p', 'picks', 'perm', 'sigma', 'noise', 'bs', 'b', 'rho', 'tol', 'x');
