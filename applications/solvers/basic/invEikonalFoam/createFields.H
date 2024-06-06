Info<< "Reading field G\n" << endl;

volScalarField G
(
    IOobject
    (
        "G",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(G)
);

// Solver error if phi is not defined.
surfaceScalarField phi
(
    IOobject(
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(G)
);

Info<< "Reading local wall thickness H from mesh\n" << endl;

volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading eikonalProperties..\n" << endl;

IOdictionary eikonalProperties
(
    IOobject
    (
        "eikonalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading reference wall thickness h from eikonalProperties\n" << endl;

dimensionedScalar h
(
    eikonalProperties.lookup("h")
);

Info<< "Reading smoothing factor sigma from eikonalProperties\n" << endl;

scalar sigma
(
    readScalar(
        eikonalProperties.lookup("sigma")
    )
);

// Variable for pseudo time term
dimensionedScalar tau
(
    eikonalProperties.lookup("tau")
);


/*
volScalarField wallDistance
(
    IOobject
    (
        "wallDistance",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("wallD", dimLength, 0.0)
);
wallDistance = wallDist(mesh).y();
*/
