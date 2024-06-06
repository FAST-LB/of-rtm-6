    if (runTime.writeTime())
    {
        volVectorField gradG(fvc::grad(G));

        volScalarField gradGx
        (
            IOobject
            (
                "gradGx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradG.component(vector::X)
        );

        volScalarField gradGy
        (
            IOobject
            (
                "gradGy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradG.component(vector::Y)
        );

        volScalarField gradGz
        (
            IOobject
            (
                "gradGz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradG.component(vector::Z)
        );

        runTime.write();
    }
