function Hbarbar = adaptivespsa_common_preconditioning(Hbar, k)
% High values of conditioner just mean that we are not using hessian
% information at all in the initial iteration when k is small.
conditioner = ((1e+2)/(k+1)^0.5);
Hbarbar = Hbar + conditioner * eye(size(Hbar));
fprintf(2, '\n conditioner %f rcond(Hbar) %f rcond(Hbarbar) %f', ...
        conditioner, rcond(Hbar), rcond(Hbarbar));