if (mesh.changing())
{	
	Info<< "Updating fiber volume fraction for new mesh." <<endl;	
    
    //const word solverTypeName = mesh.thisDb().lookupObject<motionSolver>("dynamicMeshDict").type();
	
        fiberVolFraction.ref() = fiberVolFraction.oldTime().ref()*mesh.V0()/mesh.V(); 
        forAll(fiberVolFraction.ref(),i)
        {
            if(fiberVolFraction.ref()[i]<0)
                fiberVolFraction.ref()[i]=0.0;
        }
        fiberVolFraction.correctBoundaryConditions();
       
	Info<< "Fiber Volume Fraction: max(): " << max(fiberVolFraction.internalField()).value() << " min(): " << min(fiberVolFraction.internalField()).value()<< endl;

}	
