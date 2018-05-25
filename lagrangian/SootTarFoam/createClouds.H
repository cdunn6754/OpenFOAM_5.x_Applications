Info<< "\nConstructing coal cloud" << endl;

const word cloudName = "coalCloud1";

coalCloud coalParcels
(
    cloudName,
    rho,
    U,
    g,
    slgThermo
);


// Make sure they use the analyticalSFORDevolatilization model
const word cloudFileName = cloudName + "Properties";
IOdictionary coalCloudProperties
(
     IOobject
     (
          cloudFileName,
          runTime.constant(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
     )
);

const word devolatilizationModel = 
         coalCloudProperties.subDict("subModels").lookup("devolatilisationModel");

if (devolatilizationModel != "analyticalSFORDevolatilization")
{
    FatalErrorInFunction << "You are using the " << devolatilizationModel
        << " model. Use the analyticalSFORDevolatilization model instead." 
        << " It works better with the SootTarFoam solver." << abort(FatalError);  
}

