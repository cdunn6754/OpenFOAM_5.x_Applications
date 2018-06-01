#include "readGravitationalAcceleration.H"

Info<< "Creating combustion model\n" << endl;

autoPtr<combustionModels::psiCombustionModel> combustion
(
    combustionModels::psiCombustionModel::New(mesh)
);

psiReactionThermo& baseThermo = combustion->thermo();


// Downcast thermo to fullThermo
// changes type from psiReactionThermo to sootHePsiThermo
// which enables us to use the functions in sootHePsiThermo
sootHePsiThermo<
    psiReactionThermo,
    SpecieMixture< reactingMixture< gasHThermoPhysics > >
    > &thermo = 
    dynamic_cast<
        sootHePsiThermo<
            psiReactionThermo,
            SpecieMixture<reactingMixture<gasHThermoPhysics > >
            > &> 
                                        (baseThermo);
                                       
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));

if (!composition.contains(inertSpecie))
{
    FatalErrorInFunction
        << "Specified inert specie '" << inertSpecie << "' not found in "
        << "species list. Available species:" << composition.species()
        << exit(FatalError);
}

volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);


// lagrangian effective density field - used externally (optional)
volScalarField rhoEffLagrangian
(
    IOobject
    (
        "rhoEffLagrangian",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimDensity, 0.0)
);

// dynamic pressure field - used externally (optional)
volScalarField pDyn
(
    IOobject
    (
        "pDyn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimPressure, 0.0)
);


Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

mesh.setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Set the turbulence into the combustion model
combustion->setTurbulence(turbulence());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
);


volScalarField sootQdot
(
    IOobject
    (
        "sootQdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sootQdot", dimEnergy/dimVolume/dimTime, 0.0)
);

// Create number of soot particles per kg field
volScalarField Ns
(
    IOobject
    (
        "Ns",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Create a source for the Ns equation
// this will be a source derived from the secondary tar 
// breakdown to form soot. It is already handled for the 
// Ysoot equation.
volScalarField NsSource
(
    IOobject
    (
        "NsSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("NsSource", dimensionSet(0,-3,-1, 0,0),0.0)
);

// Soot Rate Fields
volScalarField R_oxidation_O2
(
    IOobject
    (
        "R_ox_O2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_ox_O2", dimensionSet(0,-3,-1, 0,1),0.0)
);

volScalarField R_growth
(
    IOobject
    (
        "R_growth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_growth", dimensionSet(0,-3,-1, 0,1),0.0)
);

volScalarField R_oxidation_OH
(
    IOobject
    (
        "R_ox_OH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_ox_OH", dimensionSet(0,-3,-1, 0,1),0.0)
);

volScalarField R_gasification_H2O
(
    IOobject
    (
        "R_gas_H2O",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_gas_H2O", dimensionSet(0,-3,-1, 0,1),0.0)
);

volScalarField R_gasification_CO2
(
    IOobject
    (
        "R_gas_CO2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_gas_CO2", dimensionSet(0,-3,-1, 0,1),0.0)
);

volScalarField R_agglomeration
(
    IOobject
    (
        "R_agglom",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_agglom", dimensionSet(0,-3,-1, 0,1),0.0)
);

volScalarField R_kronenburg_OH
(
    IOobject
    (
        "R_Kronenburg_OH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_Kronenburg_OH", dimensionSet(0,-3,-1, 0,1),0.0)
);

// Tar Rate fields
volScalarField R_oxidation
(
    IOobject
    (
        "R_ox_tar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_ox_O2", dimensionSet(0,-3,-1, 0,1),0.0)
);

volScalarField R_cracking
(
    IOobject
    (
        "R_cr_tar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_ox_O2", dimensionSet(0,-3,-1, 0,1),0.0)
);

volScalarField R_sootFormation
(
    IOobject
    (
        "R_sf_tar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_ox_O2", dimensionSet(0,-3,-1, 0,1),0.0)
);


// Debugging field for substeps taken
volScalarField subStepsTaken
(
    IOobject
    (
        "subStep",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R_ox_O2", dimensionSet(0,0,0, 0,0),0.0)
);



// Create the system for solving the soot
TwoEquationSoot sootSystem(thermo,composition,Ns,mesh,
R_growth, R_oxidation_O2, R_oxidation_OH, 
R_gasification_H2O, R_gasification_CO2, R_agglomeration,R_kronenburg_OH);

SootTarModel tarSystem(thermo, composition, Ns, mesh,
R_oxidation,
R_cracking,
R_sootFormation,
subStepsTaken.primitiveFieldRef());

#include "createMRF.H"
#include "createClouds.H"
#include "createRadiationModel.H"
