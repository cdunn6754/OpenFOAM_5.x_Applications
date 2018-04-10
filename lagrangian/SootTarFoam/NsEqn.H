{
  
    // calculate the tar derived Ns source term
    // depending on how much soot is released from the 
    // particles
    // - 90,000 carbon atoms per soot particle
    // - 6.232e26 [atoms/kmol]
    // - 12.0 [kg/kmol] molecular weight of carbon
    NsSource.primitiveFieldRef() 
        = -coalParcels.SYi(i, Ysoot)->source() / mesh.V() * 
        (6.232e26/(9.0e4*12.01));


    fvScalarMatrix NsEqn
        (
            fvm::ddt(rho, Ns)
            + fvm::div(phi, Ns)
            - fvm::laplacian(turbulence->muEff()/SSc, Ns)
            ==
            sootSystem.sourceN()
            + NsSource
        );

    NsEqn.relax();

    fvOptions.constrain(NsEqn);

    NsEqn.solve(mesh.solver("Ns"));

    fvOptions.correct(Ns);

    // Confusing but 
    // ensures minimum value of 0.0
    Ns.max(0.0); 

}
