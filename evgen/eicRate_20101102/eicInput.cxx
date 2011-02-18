#include "eicInput.h"
#include <stdio.h>
#include <cstdlib>

eicInput::eicInput(const char *file){
    printf("Reading %s for input\n", file);

    FILE *f = fopen(file, "r");

    if( !f ){ printf("%s cannot be opened\n", file); exit(1); }

    char dummy[255];

    fscanf(f, "%s%d", dummy, &fData.nevt);
    fscanf(f, "%s%d", dummy, &fData.nprnt);
    fscanf(f, "%s%lf%s", dummy, &fData.lumin, dummy);
    fscanf(f, "%s%lf%s", dummy, &fData.runtime, dummy);
    fscanf(f, "%s%lf%s", dummy, &fData.e_energy, dummy);
    fscanf(f, "%s%lf%s", dummy, &fData.ion_energy, dummy);
    fscanf(f, "%s%lf%s", dummy, &fData.ion_mass, dummy);
    fscanf(f, "%s%d", dummy, &fData.ion_Z);
    fscanf(f, "%s%d", dummy, &fData.ion_N);
    fscanf(f, "%s%s", dummy, fData.output);

    fData.lumin *= 1e4; // Convert cm^-2 to m^-2
    fData.runtime *= 3600; // Convert hours to s

    fclose(f);
}

eicInput::~eicInput(){
}