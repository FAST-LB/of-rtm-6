        // End simulation and write last time step when cavity is filled (value can be defined in the controlDict)
	if (alpha1.weightedAverage(mesh.V()).value() > fullValue)
	{
	    Info << "I" << endl 
	    << "I" << endl
	    << "I --> End simulation" << endl 
	    << "I --> Liquid phase volume fraction = " << alpha1.weightedAverage(mesh.V()).value() << " > " << fullValue << endl
	    << "I --> Cavity Filled."<<endl;

	    runTime.writeAndEnd();
	}
	
	// End simulation and write last time step when filling gradient is too low (value can be defined in the controlDict)
	if (max(fvc::ddt(alpha1)).value() < minGradAlphaValue && runTime.time().value() > runTime.startTime().value()+runTime.deltaT().value())
	{
	    Info << "I" << endl 
	    << "I" << endl
	    << "I --> End simulation" << endl 
	    << "I --> Liquid phase filling gradient = " << max(fvc::ddt(alpha1)).value() << " < " << minGradAlphaValue << endl
	    << "I --> Cavity Filled."<<endl;

	    runTime.writeAndEnd();
	}
	
	// End simulation and write last time step when fiber volume fraction is too high (value can be defined in the controlDict)
	if (max(fiberVolFraction.internalField()).value() > maxFvfValue)
	{
            Info << "I" << endl 
	    << "I" << endl 
	    << "I --> End simulation" << endl 
	    << "I --> Maximum fiber volume fraction = " << max(fiberVolFraction.internalField()).value() << " > "<< maxFvfValue << endl
	    <<"I --> Fiber volume fraction too high."<<endl;
            
            runTime.writeAndEnd();
	}
	
if (kineticZones.active())
{
	// End simulation and write last time step when resin is cured (value can be defined in the controlDict)
	if (max(U.db().lookupObject<volScalarField>("thermo:mu").internalField()).value() > maxMuValue)
	{
	    Info << "I" << endl 
	    << "I" << endl 
	    << "I --> End simulation" << endl 
	    << "I --> Maximum viscosity = " << max(U.db().lookupObject<volScalarField>("thermo:mu").internalField()).value() << " > "<< maxMuValue << endl
	    <<"I --> Viscosity too high."<<endl;

	    runTime.writeAndEnd();
	}
}