{
    volScalarField& Yi = Y[i];

    // Tar schmidt number
    const scalar ScTar(0.7); // Taken from Brown 1998

    fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
            + mvConvection->fvmDiv(phi, Yi)
            - fvm::laplacian(turbulence->muEff()/ScTar, Yi)
            ==
            coalParcels.SYi(i, Yi)
            + combustion->R(Yi)
            + sootSystem.sourceY(Yi)
            + tarSystem.sourceY(Yi)
            + fvOptions(rho, Yi)
        );

    YiEqn.relax();

    fvOptions.constrain(YiEqn);

    YiEqn.solve(mesh.solver("Yi"));

    fvOptions.correct(Yi);

    Yi.max(0.0);
    Yt += Yi;
}
