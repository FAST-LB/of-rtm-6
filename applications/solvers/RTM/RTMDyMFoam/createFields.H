#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
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

volScalarField fiberVolFraction
(
    IOobject
    (
        "fiberVolFraction",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField KField
(
    IOobject
    (
        "permeability",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("Perm", dimLength*dimLength, Zero),
    zeroGradientFvPatchScalarField::typeName
);

volVectorField drapeDirection
(
    IOobject
    (
        "drapeDirection",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("dD", dimless, Zero)
);

volScalarField fillingTime
(
    IOobject
    (
        "fillingTime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("fT", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField resinAge
(
    IOobject
    (
        "resinAge",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rA", dimTime, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField scaleFactor
(
    IOobject
    (
        "massScaleFactor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1.0)
);

#include "createPhi.H"

Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
twoPhaseMixtureThermo mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

volScalarField mu
(
    IOobject
    (
        "muFluid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*mixture.thermo1().mu() + alpha2*mixture.thermo2().mu()
);


bool calcMu = mixture.thermo1().lookupOrDefault("calcMu", false);
if (calcMu == true)
    Info<< "Using time dependent viscosity\n" << endl;

Info<< "Reading thermophysical properties\n" << endl;

const volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& rho2 = mixture.thermo2().rho();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);


dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    mixture
);

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField dgdt(alpha1*fvc::div(phi));

#include "createAlphaFluxes.H"

// Construct compressible turbulence model
compressibleInterPhaseTransportModel turbulence
(
    rho,
    U,
    phi,
    rhoPhi,
    alphaPhi10,
    mixture
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

#include "createMRF.H"
#include "createFvOptions.H"

scalar fullValue = 0.999;
runTime.controlDict().readIfPresent("fullValue", fullValue);

scalar minGradAlphaValue = 0;
runTime.controlDict().readIfPresent("minGradAlphaValue", minGradAlphaValue);

scalar maxMuValue = 0.5;
runTime.controlDict().readIfPresent("maxMuValue", maxMuValue);

scalar maxFvfValue = 0.9;
runTime.controlDict().readIfPresent("maxFvfValue", maxFvfValue);

bool deformationMassScaling = false;
runTime.controlDict().readIfPresent("deformationMassScaling", deformationMassScaling);
