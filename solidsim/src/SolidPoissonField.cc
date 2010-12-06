#include "SolidPoissonField.hh"

SolidPoissonField::SolidPoissonField( const char *name, const char *file )
    : SolidFieldMap(name) 
{
    if( file ) fFileName = file;
    else fFileName = "";
}

/*! Add a field in file specified in file
  fFileName generated by Poisson/Superfish

  return 0 on success
  return 1 on fail
  */
int SolidPoissonField::AddField(){

    const char thisfunc[255] = "AddField";

    /// Read in, add field maps

    // OK
    return 0;

    // Some failure (file not found etc...)

    fprintf(stderr, "%s::%s Something has failed\n", GetClassName(), thisfunc);
    return 1;
} 
