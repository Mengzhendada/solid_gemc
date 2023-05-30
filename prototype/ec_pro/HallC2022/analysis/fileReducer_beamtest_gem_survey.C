#include <iostream> 
#include <fstream>
#include <cmath> 
#include "math.h" 
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TPaveText.h"
#include "TText.h"
#include "TSystem.h"
#include "TArc.h"
#include "TString.h"
#include <vector>
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TFile.h"

using namespace std;

#include "analysis_tree_solid_hgc.C"
#include "analysis_tree_solid_ec.C"
#include "analysis_tree_solid_spdsurvey.C"
#include "analysis_tree_solid_gemsurvey.C"
#define MAX_CLUSTERS_PER_PLANE 2000
#define MAX_CHANNEL 16
// some numbers to be hard coded 
// make sure they are correct while using this script
//################################################################################################################################################## 

const double DEG=180./3.1415926;   //rad to degree

//#####################################################################################################################################################

//input:
//    infile: the path of input root file from GEMC_SOLID
//    numberofFile: how many number of files in the infile, usually 1.0E4 events per file,
//                  for beam on target root tree, 1.0E9 corresponds to 80 triggers, use 80 for it
//    key:  the string used to tell what kind of run it is
//    evgen: event type, 0 is beam on target, 1 is eDIS, 2 is eAll, 3 is bggen, 4 is even file 
//int fileReducer_beamtest_gem_survey(string inputfile_name,int numberOfFile=1, double event_actual=1, int evgen=1){
int fileReducer_beamtest_gem_survey(int numberOfFile=1, double event_actual=1, int evgen=1){
 const int t=2;
        TFile *file[t];
  TTree *tree_generated[t];
  TTree *tree_flux[t];
  TTree *tree_solid_gem[t];
  TTree *tree_header[t];
  TTree *tree_solid_ec[t];
  TTree *tree_solid_ec_ps[t];
  TTree *tree_solid_hgc[t];
  TTree *tree_solid_spd[t];
   //file[0]=new TFile("/lustre19/expphy/volatile/halla/solid/tianye/container/HallC_beamtest_eAll_18deg_LD2_20m_30x30x2.5cm3_poly_beamline/HallC_beamtest_eAll_LD2_18deg_filenum9991_9.991e8_10.6GeV_shielding_30x30x2.5cm3_poly_beamline.root");
   //file[1]=new TFile("/lustre19/expphy/volatile/halla/solid/tianye/container/HallC_beamtest_eAll_18deg_LD2_20m_30x30x2.5cm3_poly_beamline/HallC_beamtest_eAll_LD2_18deg_filenum8999_8.999e8_10.6GeV_shielding_30x30x2.5cm3_poly_beamline.root");
   file[0]=new TFile("/lustre19/expphy/volatile/halla/solid/tianye/container/HallC_beamtest_bggen_LH2_pi0_18deg_LD_21m_30degcut_poly/HallC_beamtest_bggen_LH2_pi0_18deg_LD_21m_30degcut_poly_filenum42960.root");
   file[1]=new TFile("/lustre19/expphy/volatile/halla/solid/tianye/container/HallC_beamtest_bggen_LH2_pi0_18deg_LD_21m_30degcut_poly/HallC_beamtest_bggen_LH2_pi0_18deg_LD_21m_30degcut_poly_filenum6290.root");
//	char the_filename[500];
//	sprintf(the_filename, "%s",inputfile_name.substr(0,inputfile_name.rfind(".")).c_str());

//	TFile* outFile = new TFile(Form("%s_reduce_tree_rate_gem1_Etotalcut_40uA_Npe_debug.root",the_filename), "RECREATE");
	//TFile* outFile = new TFile(Form("%s_reduce_tree_rate_CC_eEcut.root",the_filename), "RECREATE");
	//TFile* outFile = new TFile("/lustre19/expphy/volatile/halla/solid/tianye/container/HallC_beamtest_eAll_18deg_LD2_20m_30x30x2.5cm3_poly_beamline/HallC_beamtest_eAll_LD2_18deg_filenum18990_10.6GeV_shielding_30x30x2.5cm3_poly_beamline_reduce_tree_rate_40uA_Npe_test.root", "RECREATE");
	TFile* outFile = new TFile("/lustre19/expphy/volatile/halla/solid/tianye/container/HallC_beamtest_bggen_LH2_pi0_18deg_LD_21m_30degcut_poly/HallC_beamtest_bggen_LH2_pi0_18deg_LD_21m_30degcut_poly_filenum49250_reduce_tree_rate_40uA_Npe_test.root", "RECREATE");
	TTree* outTree = new TTree("T", "HallC beam test simulation tree");
//	ofstream outfile1(Form("%s_output.csv",the_filename));
	cout<<"numberOfFile="<<numberOfFile<<"event_actual="<<event_actual<<"evgen="<<evgen<<endl;
	//const int ch_hgc=17;	//quad readout
	const int ch_hgc=480;	//quad readout
	const double PEthresh_hgc=1; //hgc pe threshold for each pmt
	const double PMTthresh_hgc=2; //hgc pmt threshold, at least 2pmts are fired in each sector
	Float_t npe_hgc_total=0,npe_hgc_total_trigged=0;
	float px_gen = 0;
	float py_gen = 0;
	float pz_gen = 0;
	float vx_gen = 0;
	float vy_gen = 0;
	float vz_gen = 0;
	int pid_gen = 0;

	float p_gen=0, theta_gen=0, phi_gen=0;

	float rate = 0;
	float rateRad = 0;
	float Npesum=0;
        float Cer[MAX_CHANNEL];
	int ecN=0;
	float PreShP,PreShP_e,PreShPx,PreShPy,PreShPz,PreShSum, PreSh_l, PreSh_r, PreSh_t,PreShtheta;
	//  float ShP,ShPx,ShPy,ShPz;
	float ShowerSum, Shower_l, Shower_r, Shower_t;
	float LASPD_Eendsum, SPD_Eendsum, SC_A_Eendsum, SC_D_Eendsum, SC_C_Eendsum,SC_B_Eendsum;
	float LASPD_Eend, SPD_Eend, SC_A_Eend, SC_D_Eend, SC_C_Eend,SC_B_Eend;
	float SC_A_P;
	float SC_D_P;
	float SC_C_P;
	float SC_B_P;
	float SPD_P;
	float LASPD_P;
	int GEM00_n=0,GEM10_n=0,GEM01_n=0,GEM11_n=0;
	//int Cer_n;
	int GEM00_np=0,GEM10_np=0,GEM01_np=0,GEM11_np=0;
	float GEM00_Edep,GEM10_Edep,GEM00_Etot,GEM10_Etot;
	float GEM01_Edep,GEM11_Edep,GEM01_Etot,GEM11_Etot;
	float GEM00_Edep2,GEM10_Edep2;
	float GEM00_Edep3,GEM10_Edep3;
	float GEM01_Edep2,GEM11_Edep2;
	float GEM01_Edep3,GEM11_Edep3;
	float GEM00_x[MAX_CLUSTERS_PER_PLANE],GEM00_y[MAX_CLUSTERS_PER_PLANE];
	float GEM10_x[MAX_CLUSTERS_PER_PLANE],GEM10_y[MAX_CLUSTERS_PER_PLANE];
	float GEM01_x[MAX_CLUSTERS_PER_PLANE],GEM01_y[MAX_CLUSTERS_PER_PLANE];
	float GEM11_x[MAX_CLUSTERS_PER_PLANE],GEM11_y[MAX_CLUSTERS_PER_PLANE];
	float GEM00_vx[MAX_CLUSTERS_PER_PLANE],GEM00_vy[MAX_CLUSTERS_PER_PLANE];
	float GEM10_vx[MAX_CLUSTERS_PER_PLANE],GEM10_vy[MAX_CLUSTERS_PER_PLANE];
	float GEM01_vx[MAX_CLUSTERS_PER_PLANE],GEM01_vy[MAX_CLUSTERS_PER_PLANE];
	float GEM11_vx[MAX_CLUSTERS_PER_PLANE],GEM11_vy[MAX_CLUSTERS_PER_PLANE];

	outTree->Branch("rate",       &rate,       "rate/F"      ); //vx at vertex
	outTree->Branch("vx",       &vx_gen,       "vx/F"      ); //vx at vertex
	outTree->Branch("vy",       &vy_gen,       "vy/F"      ); //vy at vertex
	outTree->Branch("vz",       &vz_gen,       "vz/F"      ); //vz at vertex
	outTree->Branch("px",       &px_gen,       "px/F"      ); //px at vertex
	outTree->Branch("py",       &py_gen,       "py/F"      ); //py at vertex
	outTree->Branch("pz",       &pz_gen,       "pz/F"      ); //pz at vertex
	outTree->Branch("p",        &p_gen,        "p/F"       ); //ptot at vertex
	outTree->Branch("pid",      &pid_gen,        "pid/I"       ); //ptot at vertex
	// preshower tree
	outTree->Branch("PreShP",      &PreShP,   "PreShP/F"    ); //momentum hit on the virtual plane of presh 
	outTree->Branch("PreShP_e",      &PreShP_e,   "PreShP_e/F"    ); //momentum hit on the virtual plane of presh 
	outTree->Branch("PreShPx",      &PreShPx,   "PreShPx/F"    ); //momentum px hit on the virtual plane of presh
	outTree->Branch("PreShPy",      &PreShPy,   "PreShPy/F"    ); //momentum py hit on the virtual plane of presh
	outTree->Branch("PreShPz",      &PreShPz,   "PreShPz/F"    ); //momentum pz hit on the virtual plane of presh
	outTree->Branch("PreShtheta",      &PreShtheta,   "PreShtheta/F"    ); //momentum pz hit on the virtual plane of presh
	outTree->Branch("PreShSum",      &PreShSum,   "PreShSum/F"    ); //Sum of the deposit energy in three preshower modules
	outTree->Branch("PreSh_l",      &PreSh_l,   "PreSh_l/F"    ); //Deposit energy in the left preshower module
	outTree->Branch("PreSh_r",      &PreSh_r,   "PreSh_r/F"    ); //Deposit energy in the right preshower module
	outTree->Branch("PreSh_t",      &PreSh_t,   "PreSh_t/F"    ); //Deposit energy in the top preshower module
	// shower tree
	outTree->Branch("ShowerSum",      &ShowerSum,   "ShowerSum/F"    ); //Sum of the deposit energy in three shower modules
	outTree->Branch("Shower_l",      &Shower_l,   "Shower_l/F"    ); //Deposit energy in the left shower module
	outTree->Branch("Shower_r",      &Shower_r,   "Shower_r/F"    ); //Deposit energy in the right shower module
	outTree->Branch("Shower_t",      &Shower_t,   "Shower_t/F"    ); //Deposit energy in the top shower module
	outTree->Branch("SC_A_P",      &SC_A_P,   "SC_A_P/F"    ); //momentum hit on the virtual plane of SC_A
	outTree->Branch("SC_A_Eendsum",      &SC_A_Eendsum,   "SC_A_Eendsum/F"    ); //total Eendosit energy in the SC_A
	outTree->Branch("SC_A_Eend",      &SC_A_Eend,   "SC_A_Eend/F"    ); //prime partile deposit energy in the SC_A
	outTree->Branch("SC_D_P",      &SC_D_P,   "SC_D_P/F"    ); //momentum hit on the virtual plane of SC_D
	outTree->Branch("SC_D_Eendsum",      &SC_D_Eendsum,   "SC_D_Eendsum/F"    ); //Deposit energy in the SC_D
	outTree->Branch("SC_D_Eend",      &SC_D_Eend,   "SC_D_Eend/F"    ); //prime particle deposit energy in the SC_D
	outTree->Branch("SC_C_P",      &SC_C_P,   "SC_C_P/F"    ); //momentum hit on the virtual plane of SC_D
	outTree->Branch("SC_C_Eendsum",      &SC_C_Eendsum,   "SC_C_Eendsum/F"    ); //Deposit energy in the SC_D
	outTree->Branch("SC_C_Eend",      &SC_C_Eend,   "SC_C_Eend/F"    ); //prime particle deposit energy in the SC_D
	outTree->Branch("SC_B_P",      &SC_B_P,   "SC_B_P/F"    ); //momentum hit on the virtual plane of SC_C
	outTree->Branch("SC_B_Eendsum",      &SC_B_Eendsum,   "SC_B_Eendsum/F"    ); //Deposit energy in the SC_C
	outTree->Branch("SC_B_Eend",      &SC_B_Eend,   "SC_B_Eend/F"    ); //prime particle deposit energy in the SC_C
	outTree->Branch("SPD_P",      &SPD_P,   "SPD_P/F"    ); //momentum hit on the virtual plane of SPD
	outTree->Branch("SPD_Eendsum",      &SPD_Eendsum,   "SPD_Eendsum/F"    ); //Deposit energy in the SPD
	outTree->Branch("SPD_Eend",      &SPD_Eend,   "SPD_Eend/F"    ); // prime particle deposit energy in the SPD
	outTree->Branch("LASPD_P",      &LASPD_P,   "LASPD_P/F"    ); //momentum hit on the virtual plane of LASPD
	outTree->Branch("LASPD_Eendsum",      &LASPD_Eendsum,   "LASPD_Eendsum/F"    ); //Deposit energy in the LASPD
	outTree->Branch("LASPD_Eend",      &LASPD_Eend,   "LASPD_Eend/F"    ); //prime partile deposit energy in the LASPD
	outTree->Branch("GEM00_n",      &GEM00_n,   "GEM00_n/I"    ); //total hits on the virtual plane of GEM00 
	outTree->Branch("GEM00_np",      &GEM00_np,   "GEM00_np/I"    ); //prime partile index hit on the virtual plane of GEM00 
	outTree->Branch("GEM00_x", &GEM00_x[0],"GEM00_x[GEM00_n]/F");// local x of GEM00
	outTree->Branch("GEM00_y", &GEM00_y[0],"GEM00_y[GEM00_n]/F");// local y of GEM00
	outTree->Branch("GEM00_vy", &GEM00_vy[0],"GEM00_vy[GEM00_n]/F");// vertex y of GEM00
	outTree->Branch("GEM00_vx", &GEM00_vx[0],"GEM00_vx[GEM00_n]/F");// vertex x of GEM00
	outTree->Branch("GEM00_Edep", &GEM00_Edep,"GEM00_Edep/F");// Deposit energy of GEM00
	outTree->Branch("GEM00_Edep2", &GEM00_Edep2,"GEM00_Edep2/F");// Deposit energy of GEM00
	outTree->Branch("GEM00_Edep3", &GEM00_Edep3,"GEM00_Edep3/F");// Deposit energy of GEM00
	outTree->Branch("GEM00_Etot", &GEM00_Etot,"GEM00_Etot/F");// Deposit energy of GEM00

	outTree->Branch("GEM10_n",      &GEM10_n,   "GEM10_n/I"    ); //total hits on the virtual plane of GEM00 
	outTree->Branch("GEM10_np",      &GEM10_np,   "GEM10_np/I"    ); //prime partile index hit on the virtual plane of GEM00 
	outTree->Branch("GEM10_x", &GEM10_x[0],"GEM10_x[GEM10_n]/F");// local x of GEM10
	outTree->Branch("GEM10_y", &GEM10_y[0],"GEM10_y[GEM10_n]/F");// local y of GEM10
	outTree->Branch("GEM10_vy", &GEM10_vy[0],"GEM10_vy[GEM10_n]/F");// vertex y of GEM10
	outTree->Branch("GEM10_vx", &GEM10_vx[0],"GEM10_vx[GEM10_n]/F");// vertex x of GEM10
	outTree->Branch("GEM10_Edep", &GEM10_Edep,"GEM10_Edep/F");// Deposit energy of GEM10
	outTree->Branch("GEM10_Edep2", &GEM10_Edep,"GEM10_Edep2/F");// Deposit energy of GEM10
	outTree->Branch("GEM10_Edep3", &GEM10_Edep,"GEM10_Edep3/F");// Deposit energy of GEM10
	outTree->Branch("GEM10_Etot", &GEM10_Etot,"GEM10_Etot/F");// Deposit energy of GEM10

	outTree->Branch("GEM01_n",      &GEM01_n,   "GEM01_n/I"    ); //total hits on the virtual plane of GEM00 
	outTree->Branch("GEM01_np",      &GEM01_np,   "GEM01_np/I"    ); //prime partile index hit on the virtual plane of GEM00 
	outTree->Branch("GEM01_x", &GEM01_x[0],"GEM01_x[GEM01_n]/F");// local x of GEM00
	outTree->Branch("GEM01_y", &GEM01_y[0],"GEM01_y[GEM01_n]/F");// local y of GEM00
	outTree->Branch("GEM01_vy", &GEM01_vy[0],"GEM01_vy[GEM01_n]/F");// vertex y of GEM00
	outTree->Branch("GEM01_vx", &GEM01_vx[0],"GEM01_vx[GEM01_n]/F");// vertex x of GEM00
	outTree->Branch("GEM01_Edep", &GEM01_Edep,"GEM01_Edep/F");// Deposit energy of GEM00
	outTree->Branch("GEM01_Edep2", &GEM01_Edep2,"GEM01_Edep2/F");// Deposit energy of GEM00
	outTree->Branch("GEM01_Edep3", &GEM01_Edep3,"GEM01_Edep3/F");// Deposit energy of GEM00
	outTree->Branch("GEM01_Etot", &GEM01_Etot,"GEM01_Etot/F");// Deposit energy of GEM00

	outTree->Branch("GEM11_n",      &GEM11_n,   "GEM11_n/I"    ); //total hits on the virtual plane of GEM00 
	outTree->Branch("GEM11_np",      &GEM11_np,   "GEM11_np/I"    ); //prime partile index hit on the virtual plane of GEM00 
	outTree->Branch("GEM11_x", &GEM11_x[0],"GEM11_x[GEM11_n]/F");// local x of GEM10
	outTree->Branch("GEM11_y", &GEM11_y[0],"GEM11_y[GEM11_n]/F");// local y of GEM10
	outTree->Branch("GEM11_vy", &GEM11_vy[0],"GEM11_vy[GEM11_n]/F");// vertex y of GEM10
	outTree->Branch("GEM11_vx", &GEM11_vx[0],"GEM11_vx[GEM11_n]/F");// vertex x of GEM10
	outTree->Branch("GEM11_Edep", &GEM11_Edep,"GEM11_Edep/F");// Deposit energy of GEM10
	outTree->Branch("GEM11_Edep2", &GEM11_Edep,"GEM11_Edep2/F");// Deposit energy of GEM10
	outTree->Branch("GEM11_Edep3", &GEM11_Edep,"GEM11_Edep3/F");// Deposit energy of GEM10
	outTree->Branch("GEM11_Etot", &GEM11_Etot,"GEM11_Etot/F");// Deposit energy of GEM10
	outTree->Branch("Npesum",      &Npesum,   "Npesum/F"    ); //Npe sum in the hgc
        outTree->Branch("Cer",&Cer[0],"Cer[16]/F");
//	outTree->Branch("Cer_n",      &Cer_n,   "Cer_n/I"    ); //total hits on the virtual plane of GEM00 
	//outTree->Branch("hit_Cer", &hit_Cer[0],"hit_Cer[Cer_n]/F");// local x of GEM10
       
	//TFile *file=new TFile(infile, "READ");
//	TFile *file=new TFile(inputfile_name.c_str());
        long int nentries[t];
	vector <double> *var1=0,*var2=0,*var3=0,*var4=0,*var5=0,*var6=0,*var7=0,*var8=0, *var9=0,*var10;
	vector <int> *gen_pid=0;
	vector <double> *gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
	vector<double> *flux_id=0,*flux_hitn=0;
	vector<double> *flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
	vector<double> *flux_trackE=0,*flux_totEdep=0,*flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0, *flux_avg_lz=0,*flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;

	TH1F *htime_photon=new TH1F("time_photon","time_photon;t (ns)",100,0,100);
  for(int n=0;n<2;n++){
	//TTree *tree_header = 0;

	tree_header[n] = (TTree*) file[n]->Get("userHeader");
	tree_header[n]->SetBranchAddress("userVar001",&var1);     //1
	tree_header[n]->SetBranchAddress("userVar002",&var2);     //x
	tree_header[n]->SetBranchAddress("userVar003",&var3);     //y
	tree_header[n]->SetBranchAddress("userVar004",&var4);     //W
	tree_header[n]->SetBranchAddress("userVar005",&var5);     //Q2
	tree_header[n]->SetBranchAddress("userVar006",&var6);     //rate
	tree_header[n]->SetBranchAddress("userVar007",&var7);     //radrate
	tree_header[n]->SetBranchAddress("userVar008",&var8);     //Ei, incoming beam energy after energy loss????
	tree_header[n]->SetBranchAddress("userVar009",&var9);     //Abeam
	tree_header[n]->SetBranchAddress("userVar010",&var10);     //Abeam

	tree_generated[n] = (TTree*) file[n]->Get("generated");
	tree_generated[n]->SetBranchAddress("pid",&gen_pid);   //particle ID
	tree_generated[n]->SetBranchAddress("px",&gen_px);     //momentum of the generated particle at target
	tree_generated[n]->SetBranchAddress("py",&gen_py);
	tree_generated[n]->SetBranchAddress("pz",&gen_pz);
	tree_generated[n]->SetBranchAddress("vx",&gen_vx);     //vertex of the generated particle at target
	tree_generated[n]->SetBranchAddress("vy",&gen_vy);
	tree_generated[n]->SetBranchAddress("vz",&gen_vz);

	tree_flux[n] = (TTree*) file[n]->Get("flux");
	tree_flux[n]->SetBranchAddress("hitn",&flux_hitn);     // hit number
	tree_flux[n]->SetBranchAddress("id",&flux_id);         //hitting detector ID
	tree_flux[n]->SetBranchAddress("pid",&flux_pid);       //pid
	tree_flux[n]->SetBranchAddress("mpid",&flux_mpid);     // mother pid
	tree_flux[n]->SetBranchAddress("tid",&flux_tid);       // track id
	tree_flux[n]->SetBranchAddress("mtid",&flux_mtid);     // mother track id
	tree_flux[n]->SetBranchAddress("otid",&flux_otid);     // original track id
	tree_flux[n]->SetBranchAddress("trackE",&flux_trackE);   //track energy of 1st step,  track here is G4 track
	tree_flux[n]->SetBranchAddress("totEdep",&flux_totEdep); //totEdep in all steps, track here is G4 track
	tree_flux[n]->SetBranchAddress("avg_x",&flux_avg_x);     //average x, weighted by energy deposition in each step
	tree_flux[n]->SetBranchAddress("avg_y",&flux_avg_y);     //average y
	tree_flux[n]->SetBranchAddress("avg_z",&flux_avg_z);     //average z
	tree_flux[n]->SetBranchAddress("avg_lx",&flux_avg_lx);   // local average x
	tree_flux[n]->SetBranchAddress("avg_ly",&flux_avg_ly);   // local average y
	tree_flux[n]->SetBranchAddress("avg_lz",&flux_avg_lz);   // local average z
	tree_flux[n]->SetBranchAddress("px",&flux_px);          // px of 1st step
	tree_flux[n]->SetBranchAddress("py",&flux_py);          // px of 1st step
	tree_flux[n]->SetBranchAddress("pz",&flux_pz);          // px of 1st step
	tree_flux[n]->SetBranchAddress("vx",&flux_vx);          // x coordinate of 1st step
	tree_flux[n]->SetBranchAddress("vy",&flux_vy);          // y coordinate of 1st step
	tree_flux[n]->SetBranchAddress("vz",&flux_vz);          // z coordinate of 1st step
	//information recorded by ec
	tree_solid_ec[n]= (TTree*) file[n]->Get("solid_ec");
	tree_solid_ec_ps[n]= (TTree*) file[n]->Get("solid_ec_ps");
	setup_tree_solid_ec(tree_solid_ec[n]);	
	setup_tree_solid_ec_ps(tree_solid_ec_ps[n]);	
	//information recorded by spd
	tree_solid_spd[n]= (TTree*) file[n]->Get("solid_spd");
	setup_tree_solid_spd(tree_solid_spd[n]);	 
	//information recorded by spd
	tree_solid_gem[n] = (TTree*) file[n]->Get("solid_gem");
	setup_tree_solid_gem(tree_solid_gem[n]);	
	tree_solid_hgc[n]= (TTree*) file[n]->Get("solid_hgc");
	setup_tree_solid_hgc(tree_solid_hgc[n]);	
	TRandom3 rand;
	rand.SetSeed(0);

	int sensor_good=0;
	int event_good=0,event_trig_good=0;
	// 	long int N_events = (long int)tree_header->GetEntries();
	nentries[n] = (long int)tree_generated[n]->GetEntries();	
        printf("Entries = %i \n",nentries[n]);

	cout << "total number of events : " << nentries[n] << endl;	

	//----------------------------
	//      loop trees
	//---------------------------
	double Ek=0;
	double Ec=0;
	double trigger_ec=0;
	int hit_id=-1,pid_max=0,mpid_max=0;
	double px_max=0,py_max=0,pz_max=0,p_max=0,p_max_e=0,theta_max=0;
	double SC_A_p_max=0;
	double SC_D_p_max=0;
	double SC_C_p_max=0;
	double SC_B_p_max=0;
	double SPD_p_max=0;
	double LASPD_p_max=0;
	int GEM00_index_max=0;
	int GEM10_index_max=0;
        int N_preshower=0;
	int pid_gen1=0;
	double hit_pf=0;
        double edep_6p1_max= 0;
		int Is_prime=-1;					  
		int Is_ECback=-1;					  
	//	double theta_gen1=0,p_gen1=0,px_gen1=0,py_gen1=0,pz_gen1=0,vx_gen1=0,vy_gen1=0,vz_gen1=0;      
	bool writeFlag = false;
	for(long int i=0;i<nentries[n];i++){	  		
		cout<<"event " << i << "\r";
		// cout<<"pid_size="<<gen_pid->size()<<endl;
		if (evgen==2){
			tree_header[n]->GetEntry(i);
			if(var5->at(0)>0.01){
			//if(var5->at(0)>1){
				//rate=(0.6667*var6->at(0))/numberOfFile; // new eDIS and eAll generator
				rate=var6->at(0)/numberOfFile; // new eDIS and eAll generator
			}
		}
		else if (evgen==0) {
			//rate=80e-6/1.6e-19/event_actual;  //beamOnTarget  file		  
			//rate=2e-6/1.6e-19/event_actual;  //beamOnTarget  file I= 5 uA	 	  
			rate=40e-6/1.6e-19/event_actual;  //beamOnTarget  file I= 5 uA	 	  
		}
		else if (evgen==4) rate=1; //even event file
		else if (evgen==3){         // bggen file
			tree_header[n]->GetEntry(i);
			rate=var10->at(0)/numberOfFile; 
		}
		else {
			cout << "Not right filemode" << endl;    
			return 0; 
		}
                //if(rate < 520833 ){
                //cout<<"event="<<i<<"rate="<<rate<<endl;
               // }
		tree_generated[n]->GetEntry(i);
		for (int j=0;j<gen_pid->size();j++) {
			pid_gen= gen_pid->at(j);
			px_gen=gen_px->at(j)/1e3;
			py_gen=gen_py->at(j)/1e3;
			pz_gen=gen_pz->at(j)/1e3;
			vx_gen=gen_vx->at(j)*0.1;
			vy_gen=gen_vy->at(j)*0.1;
			vz_gen=gen_vz->at(j)*0.1;
			p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
			//cout<<"event="<<i<<"gen_pid="<<pid_gen<<endl; 
		}
		//---
		tree_flux[n]->GetEntry(i);		  
		int sec_hgc=0;		
		int Is_trig=0;					  
		Is_prime=-1;					  
		Is_ECback=-1;					  
		hit_id=-1;
		double Eec=0,Eec_photonele=0,Eec_ele=0,EdepSC_D=0,EdepSC_C=0;
		px_max = 0;
		theta_max = 0;
		py_max = 0;
		pz_max = 0;
		p_max = 0;
		p_max_e = 0;
		pid_max = 0;
		mpid_max = 0;
		GEM00_n=0;
		GEM10_n=0;
		GEM01_n=0;
		GEM11_n=0;
                edep_6p1_max=0;
		//Cer_n=0;
		SC_A_p_max=0;
		SC_D_p_max=0;
		SC_C_p_max=0;
		SC_B_p_max=0;
		SPD_p_max=0;
		LASPD_p_max=0;
		GEM00_index_max=0;
		GEM10_index_max=0;
                N_preshower=0;
		for (Int_t j=0;j<flux_hitn->size();j++) {
			//check tid tid==1 prime particle
			//	if(flux_tid->at(j) !=1) continue;
			hit_pf=sqrt(flux_px->at(j)*flux_px->at(j)+flux_py->at(j)*flux_py->at(j)+flux_pz->at(j)*flux_pz->at(j));  //MeV to GeV
	                double hit_th=atan2(sqrt(flux_px->at(j)*flux_px->at(j)+flux_py->at(j)*flux_py->at(j)),flux_pz->at(j))*DEG;
			//double E=flux_trackE->at(j)/1e3;		  
			double E=flux_trackE->at(j);		  
			if(flux_id->at(j)==1) hit_id=0; // SC_A front 5x7 cm2
			else if(flux_id->at(j)==2) hit_id=1; // SC_D front 6.5x5.5 preshower title
			else if(flux_id->at(j)==3) hit_id=2; // EC front 
		//	else if(flux_id->at(j)==9) hit_id=9; // EC shower front
			else if(flux_id->at(j)==4) hit_id=3; // EC back 
			else if(flux_id->at(j)==5) hit_id=4; // GEM00
			else if(flux_id->at(j)==6) hit_id=5; // GEM10	  
			else if(flux_id->at(j)==7) hit_id=7; // SC4 5x10x1
			else if(flux_id->at(j)==8) hit_id=8; // LASPD	  
			else if(flux_id->at(j)==10) hit_id=6; // Cherenkov front window
			else if(flux_id->at(j)==11) hit_id=9; // SC_C 4.5x18cm2 
		//	else if(flux_id->at(j)==13) hit_id=10; // 
			else if(flux_id->at(j)==20) hit_id=10; // GEM01
			else if(flux_id->at(j)==21) hit_id=11; // GEM11
		//	else if(flux_id->at(j)==22) hit_id=13; // 
		//	else if(flux_id->at(j)==23) hit_id=14; // 
		//	else if(flux_id->at(j)==24) hit_id=15; // 
		//	else if(flux_id->at(j)==30) hit_id=16; //
			else cout << "wrong flux_id" << flux_id->at(j) << endl;
				if(hit_id==4 /*&& flux_pid->at(j)==22 && flux_tid->at(j)==1 && flux_mtid->at(j)==0*/){
                                Is_prime=1;
		                N_preshower ++;
				if (E >= edep_6p1_max){
					edep_6p1_max=E;
					p_max=hit_pf;
					px_max=flux_px->at(j);
					py_max=flux_py->at(j);
					pz_max=flux_pz->at(j);
                                        theta_max = hit_th;
                                        pid_max = flux_pid->at(j);
                                        mpid_max = flux_mpid->at(j);
                                 }
				} 
				if(hit_id==0 ){
					SC_A_p_max=hit_pf;
				}
				if(hit_id==1){
					SC_D_p_max=hit_pf;
				}
				if(hit_id==7){
                                Is_ECback=1;
					SC_B_p_max=hit_pf;
				}
				if(hit_id==9){
					SC_C_p_max=hit_pf;
				}
				//if(hit_id==7){
				//	SPD_p_max=hit_pf;
				//}
				if(hit_id==8){
					LASPD_p_max=hit_pf;
				}
				if(hit_id==4){
					GEM00_np=j;
				}
				if(hit_id==5){
					GEM10_np=j;
				}
				if(hit_id==10){
					GEM01_np=j;
				}
				if(hit_id==11){
					GEM11_np=j;
				}
			//}// end tid==1
			if(hit_id==4){
				GEM00_x[GEM00_n]=flux_avg_lx->at(j)*0.1;
				GEM00_y[GEM00_n]=flux_avg_ly->at(j)*0.1;
				GEM00_vx[GEM00_n]=flux_avg_x->at(j)*0.1;
				GEM00_vy[GEM00_n]=flux_avg_y->at(j)*0.1;
				GEM00_n++;
			}
			if(hit_id==5){
				GEM10_x[GEM10_n]=flux_avg_lx->at(j)*0.1;
				GEM10_y[GEM10_n]=flux_avg_ly->at(j)*0.1;
				GEM10_vx[GEM10_n]=flux_avg_x->at(j)*0.1;
				GEM10_vy[GEM10_n]=flux_avg_y->at(j)*0.1;
				GEM10_n++;
			}
			if(hit_id==10){
				GEM01_x[GEM01_n]=flux_avg_lx->at(j)*0.1;
				GEM01_y[GEM01_n]=flux_avg_ly->at(j)*0.1;
				GEM01_vx[GEM01_n]=flux_avg_x->at(j)*0.1;
				GEM01_vy[GEM01_n]=flux_avg_y->at(j)*0.1;
				GEM01_n++;
			}
			if(hit_id==11){
				GEM11_x[GEM11_n]=flux_avg_lx->at(j)*0.1;
				GEM11_y[GEM11_n]=flux_avg_ly->at(j)*0.1;
				GEM11_vx[GEM11_n]=flux_avg_x->at(j)*0.1;
				GEM11_vy[GEM11_n]=flux_avg_y->at(j)*0.1;
				GEM11_n++;
			}
			if(hit_id==6 && abs(flux_avg_lx->at(j)*0.1)<2 && abs(flux_avg_ly->at(j)*0.1)<2){                                
				Is_trig=1;
			}
			//cout<<"event="<<i<<"j="<<j<<" GEM00_n="<<GEM00_n<<"GEM10_n="<<GEM10_n<<endl;
		}	// end of flux		
                //cout<<"event="<<i<<"preshower_N="<<N_preshower<<endl;
                // process gem
           //   if(Is_prime==1 &&  N_preshower<=1 && edep_6p1_max==200){

		tree_solid_gem[n]->GetEntry(i);
		double totEdep_gem1_gas[6]={0};
		double totEdep_gem2_gas[6]={0};
		double totEdep_gem1tot=0;
		double totEdep_gem2tot=0;
		double totEdep_gem3_gas[6]={0};
		double totEdep_gem4_gas[6]={0};
		double totEdep_gem3tot=0;
		double totEdep_gem4tot=0;
		int procID_gem=0;
		int mpid_gem=0;
		int pid_gem=0;
	//		process_tree_solid_gem(tree_solid_gem,totEdep_gem1_gas,totEdep_gem2_gas,totEdep_gem1tot,totEdep_gem2tot, procID_gem,mpid_gem,pid_gem);
			process_tree_solid_gem(tree_solid_gem[n],totEdep_gem1_gas,totEdep_gem2_gas,totEdep_gem3_gas,totEdep_gem4_gas,totEdep_gem1tot,totEdep_gem2tot,totEdep_gem3tot,totEdep_gem4tot, procID_gem,mpid_gem,pid_gem);
                GEM00_Edep=totEdep_gem1_gas[1];
                GEM01_Edep=totEdep_gem2_gas[1];
                GEM10_Edep=totEdep_gem3_gas[1];
                GEM11_Edep=totEdep_gem4_gas[1];
                GEM00_Etot=totEdep_gem1tot;
                GEM10_Etot=totEdep_gem2tot;
                GEM01_Etot=totEdep_gem3tot;
                GEM11_Etot=totEdep_gem4tot;
                GEM00_Edep2=totEdep_gem1_gas[2];
                GEM10_Edep2=totEdep_gem2_gas[2];
                GEM00_Edep3=totEdep_gem1_gas[3];
                GEM10_Edep3=totEdep_gem2_gas[3];
                GEM01_Edep2=totEdep_gem3_gas[2];
                GEM11_Edep2=totEdep_gem4_gas[2];
                GEM01_Edep3=totEdep_gem3_gas[3];
                GEM11_Edep3=totEdep_gem4_gas[3];

		// process ec
		tree_solid_ec[n]->GetEntry(i);
		tree_solid_ec_ps[n]->GetEntry(i);
		//if(hit_preshowerN>0){
		double Eend_ec_sum=0;
		double Eend_ec_ps_sum=0;	
		double Eend_ec[4]={0};
		double Eend_ec_ps[4]={0};
		process_tree_solid_ec(tree_solid_ec[n],tree_solid_ec_ps[n],Eend_ec_sum,Eend_ec_ps_sum,Eend_ec,Eend_ec_ps);
		PreShSum=Eend_ec_ps_sum;
		PreSh_l=Eend_ec_ps[2];
		PreSh_r=Eend_ec_ps[3];
		PreSh_t=Eend_ec_ps[1];
		ShowerSum=Eend_ec_sum;
		Shower_l=Eend_ec[2];
		Shower_r=Eend_ec[3];
		Shower_t=Eend_ec[1];
		//PreShP=p_max;
		PreShP=edep_6p1_max;
		PreShPx=px_max;
		PreShPy=py_max;
		PreShPz=pz_max;
		PreShtheta=theta_max;
		tree_solid_spd[n]->GetEntry(i);
		double Edep_sc1=0, Edep_sc2=0, Edep_sc3=0, Edep_sc4=0, Edep_spd1=0, Edep_spd2=0;
		process_tree_solid_spd_simple(tree_solid_spd[n],Edep_sc1,Edep_sc2,Edep_sc3,Edep_sc4,Edep_spd1,Edep_spd2);
		//SC_A
		SC_A_P=SC_A_p_max;
		SC_A_Eendsum = Edep_sc1;
		if(SC_A_P>0){
			SC_A_Eend = Edep_sc1;
			//cout<<"event="<<i<<"SC_A_P="<<SC_A_P<<"SC_A_Eendsum="<<SC_A_Eendsum<<endl;
		}
		//SC_D
		SC_D_P=SC_D_p_max;
		SC_D_Eendsum = Edep_sc2;
		if(SC_D_P>0){
			SC_D_Eend = Edep_sc2;
		}
                //SC_C
		SC_C_P=SC_C_p_max;
		SC_C_Eendsum = Edep_sc3;
		if(SC_C_P>0){
			SC_C_Eend = Edep_sc3;
		}
                //SC_B
		SC_B_P=SC_B_p_max;
		SC_B_Eendsum = Edep_sc4;
		if(SC_B_P>0){
			SC_B_Eend = Edep_sc4;
		}
		//SPD
		SPD_P=SPD_p_max;
		SPD_Eendsum = Edep_spd1;
		if(SPD_Eendsum>0){
			SPD_Eend = Edep_spd1;
		}
		//LASPD
		LASPD_P=LASPD_p_max;
		LASPD_Eendsum = Edep_spd2;
		if(LASPD_Eendsum>0){
			LASPD_Eend = Edep_spd2;
		}
		//--- hgc 
		//-----------------------
		tree_solid_hgc[n]->GetEntry(i);
		double hit_hgc[ch_hgc]={0};
		int trigger_hgc[30]={0};
		int ntrigsecs_hgc=0;
		int ntrigpmts_hgc=0;
		int photon_mtid=0;
                process_tree_solid_hgc(tree_solid_hgc[n],hit_hgc,trigger_hgc,ntrigsecs_hgc,PMTthresh_hgc,PEthresh_hgc,ch_hgc,htime_photon,photon_mtid);
		npe_hgc_total=0;
		npe_hgc_total_trigged=0;		  
		Cer[ch_hgc]={0};
		for(int sensor_id=0;sensor_id<ch_hgc;sensor_id++) {	
                                Cer[sensor_id] = hit_hgc[sensor_id];  
			if (hit_hgc[sensor_id]>0){
				npe_hgc_total += hit_hgc[sensor_id]; 
                                //hit_Cer[sensor_id]=hit_hgc[sensor_id];  
                        //cout<<"event="<<i<<"sensor_id="<<sensor_id<<"Cer="<<hit_Cer[sensor_id]<<"hit_hgc="<<hit_hgc[sensor_id]<<"CerSum="<<npe_hgc_total<<endl;
			}
		}
              //if(Is_trig==1){
              //if(npe_hgc_total>0){
               // cout<<"event="<<i<<"CerSum="<<npe_hgc_total<<endl;
              // }
		Npesum=npe_hgc_total;
		outTree->Fill();
          //}//prime cut
	} //end loop
	//do outputs
      }//file loop
	cout<<endl;
	outFile->cd();
	outTree->Write();
	outFile->Close();
	}
