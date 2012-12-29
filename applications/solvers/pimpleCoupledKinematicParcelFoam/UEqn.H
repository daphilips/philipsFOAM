// Solve the momentum equation
fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
  ==
    sources(U) 
  + irhoVal*kinematicCloud.SU(U)
);

UEqn.relax();

sources.constrain(UEqn);


if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
        ==
        fvc::reconstruct
        (
          -fvc::snGrad(p_rgh)*mesh.magSf()
        )
    );
}
