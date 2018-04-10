{

    // Soot mass fraction
    volScalarField& Ysoot = Y[i];

    // Soot schimdt number
    const scalar SSc(700.0); // from BYU paper

    fvScalarMatrix YsootEqn
        (
            fvm::ddt(rho, Ysoot)
            + mvConvection->fvmDiv(phi, Ysoot)
            - fvm::laplacian(turbulence->muEff()/SSc, Ysoot)
            ==
            sootSystem.sourceY(Ysoot)
            + coalParcels.SYi(i, Ysoot)
            + combustion->R(Ysoot)
            + fvOptions(rho, Ysoot)
        );
    
    #include "NsEqn.H"
   
    YsootEqn.relax();

    fvOptions.constrain(YsootEqn);

    YsootEqn.solve(mesh.solver("Yi"));

    fvOptions.correct(Ysoot);

    Ysoot.max(0.0);
    Yt += Ysoot;
}
