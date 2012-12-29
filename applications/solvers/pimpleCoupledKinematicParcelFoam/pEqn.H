{
  volScalarField rAU("rAU", 1.0/UEqn.A());
  surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));

  U = rAU*UEqn.H();

  phi = (fvc::interpolate(U) & mesh.Sf())
       + fvc::ddtPhiCorr(rAU, U, phi);

  adjustPhi(phi, U, p_rgh);

  while (pimple.correctNonOrthogonal())
  {
      fvScalarMatrix p_rghEqn
      (
          fvm::laplacian(rAU, p_rgh) == fvc::div(phi)
      );

      p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

      p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

      if (pimple.finalNonOrthogonalIter())
      {
          // Calculate the conservative fluxes
          phi -= p_rghEqn.flux();

          // Explicitly relax pressure for momentum corrector
          p_rgh.relax();

          // Correct the momentum source with the pressure gradient flux
          // calculated from the relaxed pressure
          U -= rAU*fvc::reconstruct((p_rghEqn.flux())/rAUf);
          U.correctBoundaryConditions();
      }
  }

  #include "continuityErrs.H"

  p = rhoVal*(p_rgh + gh);

  if (p_rgh.needReference())
  {
      p += dimensionedScalar
      (
          "p",
          p.dimensions(),
          pRefValue - getRefCellValue(p, pRefCell)
      );
      p_rgh = irhoVal*p - gh;
  }
}
