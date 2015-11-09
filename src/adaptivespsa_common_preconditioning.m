function Hbarbar = adaptivespsa_common_preconditioning(Hbar, k)
delta_eye = 1e-16 * exp(-k) * eye(size(Hbar));
Hbarbar = sqrtm(Hbar * Hbar + delta_eye);
% Hbarbar = Hbar + delta_eye; % This conditioning is incredibly bad!
fprintf(2, '\n rcond(Hbar) %f rcond(Hbarbar) %f', ...
        rcond(Hbar), rcond(Hbarbar));