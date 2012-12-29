    Info<< "\nReading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    dimensionedScalar rhoVal
    (
        transportProperties.lookup("rho")
    );
    Info<< "rhoVal is " << rhoVal << endl;

    dimensionedScalar irhoVal("irhoVal",(1.0/rhoVal));
