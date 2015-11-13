function Bbar = rank_two_update_v2(Bbar, a, b, delta_tilda_k, delta_k)
% Stage 1 update
tmp_stage = Bbar * delta_tilda_k;
tmp_deno = a + b * (delta_k' * tmp_stage);
Bbar = Bbar - tmp_stage * ((delta_k' *(b/tmp_deno)) * Bbar);
% Stage 2 update
tmp_stage = Bbar * delta_k;
tmp_deno = a + b * (delta_tilda_k' * tmp_stage);
Bbar = Bbar - tmp_stage * ((delta_tilda_k'*(b/tmp_deno)) * Bbar);
