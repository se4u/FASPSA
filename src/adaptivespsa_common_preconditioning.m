function Hbarbar = adaptivespsa_common_preconditioning(Hbar, k)
delta_eye = 1e-8 * exp(-k) * eye(size(Hbar));
Hbarbar = sqrtm(Hbar * Hbar + delta_eye);
fprintf(2, '\n rcond(Hbar) %f rcond(Hbarbar) %f', ...
        rcond(Hbar), rcond(Hbarbar));