    Info<< "Reading normalised dynamic presure field p_rgh\n" << endl;
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

    Info<< "Creating phi [m3/s] " << endl;
    #include "createPhi.H"

    Info<< "Creating turbulence model\n" << endl;
    singlePhaseTransportModel laminarTransport(U, phi);
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    Info<< "Creating source list\n" << endl;
    IObasicSourceList sources(mesh);


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    Info<< "Creating absolute pressure field p [Pa]\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rhoVal*(p_rgh + gh)
    );

    Info<< "Reading field U\n" << endl;

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }

    Info<< "Creating mass flux density surface field rhoPhi [kg/s] " << endl;
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        rhoVal*linearInterpolate(U) & mesh.Sf()
    );

    // dummy-fields to satisfy the requirements of the kinematicCloud
    Info<< "Creating rho field " << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        rhoVal,
        "zeroGradient"
    );

    Info<< "Creating mu field " << endl;
    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        rhoVal*nu,
        "zeroGradient"
    );

    Info<< "Constructing kinematicCloud" << endl;
    basicKinematicCollidingCloud kinematicCloud
    (
        "kinematicCloud",
        rho,
        U,
        mu,
        g
    );

