sgm = alpha1 * sgm1 + (scalar(1) - alpha1) * sgm2;
eps = scalar(1) / (alpha1 / eps1 + (scalar(1) - alpha1) / eps2);

surfaceScalarField rhoEFlux (-fvc::interpolate(sgm)*mesh.magSf()*fvc::snGrad(Ue));

fvScalarMatrix rhoEEqn
(
    fvm::ddt(rhoE)
  + fvm::div(phi, rhoE)
  + fvc::div(rhoEFlux)
);
rhoEEqn.solve();

// solve
fvScalarMatrix UeEqn
(
    fvm::laplacian(eps,Ue) == -rhoE
);

UeEqn.solve();
