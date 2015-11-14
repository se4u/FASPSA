clc; clear;
p = 3;
% Set up the temporary storage once.
time_old = 0;
time_new = 0;
for i = 1:1
    Bbar = eye(p);
    Bbar = (Bbar + Bbar');
    b = 2;
    a = 3;
    delta_tilda_k = rand(p, 1);
    delta_k = rand(p, 1);

    tic
    Bbar_gold = rank_two_update_v2(Bbar, a, b, delta_tilda_k, delta_k);
    time_old = time_old + toc;

    tmp_stage = zeros(p, 1);
    Binv = zeros(p);
    tmp_stage_b = zeros(p, 1);


    tic
    Bbar_fast = rank_two_update_v2_fast(...
        Bbar, a, b, delta_tilda_k, delta_k, ...
        tmp_stage, Binv, tmp_stage_b);
    time_new = time_new + toc;

    % (all(all(Bbar_fast == Bbar_gold)))
    norm(Bbar_fast - Bbar_gold)
end
[time_old time_new]