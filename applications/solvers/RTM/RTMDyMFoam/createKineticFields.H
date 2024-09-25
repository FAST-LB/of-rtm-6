Info<< "Reading field cure\n" << endl;
volScalarField& cure = mixture.cure();
Info<< "Reading field cureRate\n" << endl;
volScalarField cureRate
(
    IOobject
    (
        "cureRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar("zero", dimless/dimTime, scalar(0.0))
);
Info<< "Reading field cureEff\n" << endl;

volScalarField cureEff
(
    IOobject
    (
        "cureEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar("zero", dimless, scalar(0.0))
);

IOkineticModelList kineticZones(mesh);

Info<< "Reading field Tg\n" << endl;
volScalarField Tg
(
    IOobject
    (
        "Tg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
	dimensionedScalar("zero", dimTemperature, scalar(0.0))
);

volScalarField TgEff
(
    IOobject
    (
        "TgEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar("zero", dimTemperature, scalar(0.0))
);

IOTgModelList TgZones(mesh);

// material state (0=liquid,1=rubbery,2=glassy)
volScalarField materialState
(
    IOobject
    (
        "materialState",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
	dimensionedScalar("zero", dimless, scalar(0.0))
);
