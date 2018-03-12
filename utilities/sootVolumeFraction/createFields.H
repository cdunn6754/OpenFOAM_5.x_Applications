// thermo model. Problably should figure out how to do this without
// making the whole combustion model
autoPtr<combustionModels::psiCombustionModel> combustion
(
    combustionModels::psiCombustionModel::New(mesh)
);
psiReactionThermo& thermo = combustion->thermo();

// Create number of soot particles per kg field
volScalarField Ns
(
    IOobject
    (
        "Ns",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// Read the soot mass fraction field
volScalarField Ysoot
(
    IOobject
    (
        "SOOT",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// Read the soot mass fraction field
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


// Field to hold the particle diameters
volScalarField Dp
(
    IOobject
    (
        "Dp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// Field to hold the soot volume fraction
volScalarField Vsoot
(
    IOobject
    (
        "Vsoot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);
