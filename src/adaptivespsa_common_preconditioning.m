function Hbarbar = adaptivespsa_common_preconditioning(Hbar, k)
Hbarbar = Hbar + ((1e-4)/(k+1)) * speye(size(Hbar));