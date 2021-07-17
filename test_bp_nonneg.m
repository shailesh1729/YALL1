clear; clc;
rng('default');

% problem sizes 
n = 1000; m = 300; k = 60; 
A = randn(m,n);
xs = zeros(n,1);
p = randperm(n);
xs(p(1:k)) = randn(k,1);
% make sure that entries in x are all non-negative
xs = abs(xs);
b = A*xs;
[Q, R] = qr(A',0);
A = Q'; b = R'\b;
opts.tol = 5e-3; 
opts.delta = 0;
opts.rho = 0;
opts.nonneg = 1; 
opts.print = 1; 

t0 = tic; 
[x,Out] = yall1(A, b, opts);
err = norm(x-xs)/norm(xs);
T = toc(t0);
fprintf('\nrelative error = %e, time %6.2f ms\n',err, T*1000);
disp(Out)
save(sprintf('data/partial_onb_bp_nonneg.mat'), 'A', 'b', 'x', 'xs');
