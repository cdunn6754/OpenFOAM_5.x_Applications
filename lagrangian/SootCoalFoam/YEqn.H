tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    Info <<"Correcting combustion " << endl;
    combustion->correct();
    sootSystem.updateSources();
    sootQdot = sootSystem.Qdot();
    Qdot = combustion->Qdot();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if ((Y[i].name() != inertSpecie) && (Y[i].name() !="SOOT"))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
                (
                    fvm::ddt(rho, Yi)
                    + mvConvection->fvmDiv(phi, Yi)
                    - fvm::laplacian(turbulence->alphaEff(), Yi)
                    ==
                    coalParcels.SYi(i, Yi)
                    + combustion->R(Yi)
                    + sootSystem.sourceY(Yi)
                    + fvOptions(rho, Yi)
                );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
      else if (Y[i].name() == inertSpecie)
        {
            inertIndex = i;
        }
      else if (Y[i].name() == "SOOT")
	{
            #include "Ysoot.H"
	}
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
