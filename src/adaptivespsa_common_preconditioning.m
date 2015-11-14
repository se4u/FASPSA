function Hbarbar = adaptivespsa_common_preconditioning(Hbar, k)
delta_eye = 1e-8;
Hbar_sq = Hbar * Hbar;
Hbarbar_f  = @(regularizer) sqrtm(Hbar_sq + regularizer * eye(size(Hbar)));
Hbarbar = Hbarbar_f(delta_eye);
while not(all(isreal(Hbarbar)))
    % exit(1);
    delta_eye = 10 * delta_eye;
    my_fprintf(2, '\n delta_eye %g ', delta_eye);
    Hbarbar = Hbarbar_f(delta_eye);
end