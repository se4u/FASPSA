clc; clear;
mex('-v', '-largeArrayDims', '-lmwblas', '-DSAFE', 'rank_two_update_v2_fast.c');
rand('seed',31415927);
randn('seed',3111113);
p_arr = 10:10:200;
p_arr = 10;
time_old_arr = [];
time_new_arr = [];
ndiff = 0;
nsym_old = 0;
nsym_new = 0;
for p = p_arr
    time_old = 0;
    time_new = 0;
    trials = 100;
    for i = 1:trials
        Bbar = rand(p);
        Bbar = (Bbar + Bbar');
        b = 2;
        a = 3;
        delta_tilda_k = rand(p, 1);
        delta_k = rand(p, 1);

        tic
        Bbar_gold = rank_two_update_v2(Bbar, a, b, delta_tilda_k, delta_k);
        time_old = time_old + toc;
        nsym_old = max(nsym_old, max(max(abs(Bbar_gold - Bbar_gold'))));

        tic
        Bbar_fast = rank_two_update_v2_fast(...
            Bbar, a, b, delta_tilda_k, delta_k);
        time_new = time_new + toc;
        nsym_new = max(nsym_new, max(max(abs(Bbar_fast - Bbar_fast'))));
        ndiff = max(ndiff, norm(Bbar_fast - Bbar_gold));
    end
    time_new_arr = [time_new_arr time_new/trials];
    time_old_arr = [time_old_arr time_old/trials];
end
fprintf(2, ...
        '\n ndiff %g, nsym_old %g, nsym_new %g ', ...
        ndiff, nsym_old, nsym_new);
if ~usejava('jvm')
    exit(~(ndiff < 1e-10));
else
    figure();
    plot(p_arr, time_old_arr, 'b');
    hold on;
    plot(p_arr, time_new_arr, 'r');
    legend('time_old_arr', 'time_new_arr');
end