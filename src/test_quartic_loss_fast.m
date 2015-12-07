p = 200; l = quartic_loss_factory(p);
time_slow = 0;
time_fast = 0;
runs = 1000;
slow = zeros(runs, 1);
fast = zeros(runs, 1);
a = randn(p, runs);
tic
for idx=1:runs
    slow(idx) = l(a(:, idx));
end
time_slow = time_slow + toc;

tic
for idx=1:runs
    fast(idx) = quartic_loss_fast(a(:, idx));
end
time_fast = time_fast + toc;

assert(max(abs(slow - fast)) < 1e-13);
fprintf(2, '\n time_slow/time_fast %e ', time_slow/time_fast);
