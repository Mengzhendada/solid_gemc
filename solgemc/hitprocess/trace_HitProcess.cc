// %%%%%%%%%%
// G4 headers
// %%%%%%%%%%
#include "G4UnitsTable.hh"


// %%%%%%%%%%%%%
// gemc headers
// %%%%%%%%%%%%%
#include "trace_HitProcess.h"
#include "SolPrimaryGeneratorAction.h"

trace_HitProcess :: trace_HitProcess() {
    // PrimaryGeneratorAction can contain more information
    // that we want to pull out, such as info on the primary
    // track, any weighting that we want
    fPMG = NULL;
}

SolPrimaryGeneratorAction *trace_HitProcess::GetPMG(){
    fPMG = SolPrimaryGeneratorAction::GetInstance();
    return fPMG;
}

PH_output trace_HitProcess :: ProcessHit(MHit* aHit, gemc_opts Opt){
	PH_output out;
	out.identity = aHit->GetId();
	HCname = "Trace Hit Process";

	if( !fPMG ){GetPMG();}

	double weight;
	if( fPMG ){ 
	    weight = fPMG->GetWeight(); 
	} else {
	    weight = 1.0;
	}

	// %%%%%%%%%%%%%%%%%%%
	// Raw hit information
	// %%%%%%%%%%%%%%%%%%%
	int nsteps = aHit->GetPos().size();

	// Get Total Energy deposited
	double Etot = 0;
	vector<G4double> Edep = aHit->GetEdep();
	for(int s=0; s<nsteps; s++){
	    Etot = Etot + Edep[s];
	}
	
	// average global, local positions of the hit
	double x, y, z;
	double lx, ly, lz;
	x = y = z = lx = ly = lz = 0;
	vector<G4ThreeVector> pos  = aHit->GetPos();
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	G4ThreeVector p = aHit->GetMom();
	
	for(int s=0; s<nsteps; s++) {
	    x  = x  +  pos[s].x()/((double) nsteps);
	    y  = y  +  pos[s].y()/((double) nsteps);
	    z  = z  +  pos[s].z()/((double) nsteps);
	    lx = lx + Lpos[s].x()/((double) nsteps);
	    ly = ly + Lpos[s].y()/((double) nsteps);
	    lz = lz + Lpos[s].z()/((double) nsteps);
	}
	
	// average time
	double time = 0;
	vector<G4double> times = aHit->GetTime();
	for(int s=0; s<nsteps; s++) time = time + times[s]/nsteps;

	// Energy of the track
	double Ene = aHit->GetE();

	out.raws.push_back(Etot/MeV);
	out.raws.push_back(x/mm);
	out.raws.push_back(y/mm);
	out.raws.push_back(z/mm);
	out.raws.push_back(lx/mm);
	out.raws.push_back(ly/mm);
	out.raws.push_back(lz/mm);
	out.raws.push_back(time/ns);
	out.raws.push_back((double) aHit->GetPID());
	out.raws.push_back(aHit->GetVert().getX()/mm);
	out.raws.push_back(aHit->GetVert().getY()/mm);
	out.raws.push_back(aHit->GetVert().getZ()/mm);
	out.raws.push_back(Ene/MeV);
	out.raws.push_back((double) aHit->GetmPID());
	out.raws.push_back(aHit->GetmVert().getX()/mm);
	out.raws.push_back(aHit->GetmVert().getY()/mm);
	out.raws.push_back(aHit->GetmVert().getZ()/mm);
	out.raws.push_back(aHit->GetTId());
	out.raws.push_back(weight);
	out.raws.push_back(p.x()/MeV);
	out.raws.push_back(p.y()/MeV);
	out.raws.push_back(p.z()/MeV);

	int id  = out.identity[0].id;
	out.dgtz.push_back(id);

	return out;
}

vector<identifier>  trace_HitProcess :: ProcessID(vector<identifier> id, G4Step* aStep, detector Detector, gemc_opts Opt)
{
    return id;
}











