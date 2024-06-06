forAll(mesh.C(), cellI)
{
    if ( alpha1[cellI] < 0.5 )
    {
        fillingTime[cellI] = runTime.time().value();        
    }    
}

forAll(mesh.C(), cellI)
{
    if ( alpha1[cellI] < 0.5 )
    {
        resinAge[cellI] += runTime.deltaT().value();        
    }    
}
