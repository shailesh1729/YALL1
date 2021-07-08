clear; clc;
rng('default');

% problem sizes 
n = 1000; m = 300; k = 60; 
sigma = 0.01;

A = randn(m,n);
xs = zeros(n,1);
p = randperm(n);
xs(p(1:k)) = randn(k,1);
noise = sigma*randn(m,1);
bs = A*xs;
b = bs + noise;
delta = norm(noise);
fprintf('sigma = %6.3e, delta = %6.3e\n\n',sigma,delta);
[Q, R] = qr(A',0);
A = Q'; b = R'\b;
opts.tol = 1e-3; 
opts.delta = delta;
opts.print = 1; 

t0 = tic; 
[x,Out] = yall1(A, b, opts);
err = norm(x-xs)/norm(xs);
T = toc(t0);
fprintf('\nrelative error = %e, time %6.2f ms\n',err, T*1000);
disp(Out)
save(sprintf('data/partial_onb_l1_l2con.mat'), 'A', 'b', 'bs', 'x', 'xs');
