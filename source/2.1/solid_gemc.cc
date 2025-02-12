/// \mainpage
/// \htmlonly <center><img src="gemc_logo.gif" width="130"></center>\endhtmlonly
/// \section overview Overview
/// gemc (<b>GE</b>ant4 <b>M</b>onte<b>C</b>arlo) GEMC is a C++ framework
/// based on <a href="http://geant4.web.cern.ch/geant4/"> Geant4 </a>
/// Libraries to simulate the passage of particles through matter.\n
/// The simulation parameters are external to the software:
/// Geometry, Materials, Fields, Banks definitions are stored in
/// external databases in various format and are loaded at run
/// time using factory methods.\n
/// \section databases Databases
/// gemc currently supports <i> mysql, gdml, TEXT </i> to build the detector systems: \n
/// - Geometry.
/// - Sensitive Detectors.
/// - Hit Process.
/// - Banks Format.
/// - Materials.
/// \section platforms Platforms Supported:
/// - <i> Windows 7 to come October 2013 </i>
/// - Linux (32, 64)
/// - Mac OS X
/// \section docs Documentation:
/// - <a href="http://gemc.jlab.org">  gemc website </a>
/// \section Dependencies:
/// - geant4 (simulation libraries)
/// - clhep  (random generators, physics vectors, geometry and linear algebra libraries)
/// - qt4 (graphic libraries)
/// - mysql
/// - scons (build system)
/// - xercesc
/// \image html gemc_logo.gif
/// \n\n
/// \author \n &copy; Maurizio Ungaro
/// \author e-mail: ungaro@jlab.org\n\n\n
/// \file gemc.cc
/// Defines the gemc main( int argc, char **argv )
/// \author \n &copy; Maurizio Ungaro
/// \author e-mail: ungaro@jlab.org\n\n\n

const char *GEMC_VERSION = "gemc 2.1";

// G4 headers
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4VisExecutive.hh"
#include "G4VModularPhysicsList.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"
#include "G4UIQt.hh"
#include "G4Qt.hh"

// Qt headers
#include <QApplication>
#include <QSplashScreen>

// gemc headers
#include "HitProcess_MapRegister.h"
#include "detector_factory.h"
#include "gemc_MainGui.h"
#include "gbank.h"
#include "MDetectorConstruction.h"
#include "MEventAction.h"
#include "outputFactory.h"
#include "HitProcess.h"
#include "PhysicsList.h"
#include "MPrimaryGeneratorAction.h"
#include "MSteppingAction.h"
#include "options.h"
#include "dmesg_init.h"
#include "run_conditions.h"
#include "fieldFactory.h"
#include "material_factory.h"
#include "mirrors_factory.h"
#include "parameter_factory.h"
#include "string_utilities.h"
#include "utils.h"

// c++ headers
#include <unistd.h>  // needed for get_pid

//header for addtional hitprocess
#include "solid_hitprocess.h"

/////////////////////////
/// <b> Main Program </b>
/////////////////////////
///  -# Sets the goptions\n
///  -# Starts QT application if USE_GUI=1
///  -# Starts the CLHEP random engine
///  -# Instantiates the Geant4 Run Manager
///  -# Builds detector map object from database
///  -# Builds Processes Routines Map
///  -# Builds Materials Map
///  -# Builds G4 Physical Volumes
///  -# Initialize Physics List
///  -# Initialize Generator
///  -# Initialize Event Action
///  -# Initialize G4Qt User Interface if USE_GUI>0
///  -# Initialize Visualization Manager if USE_GUI>0


// get_pid is useful only on the farm to set the seed
// can set to zero in Windows environment
// ideally we'd want __get_pid();
#ifdef _MSC_VER
#include <stdio.h>
#include <process.h>
	//int get_pid(){return __get_pid();}
	int get_pid(){return 0;}
#endif



int main( int argc, char **argv )
{
	clock_t startTime = clock();
	cout << endl;
	goptions gemcOpt;
	gemcOpt.setGoptions();
	gemcOpt.setOptMap(argc, argv);
	
	double use_gui   = gemcOpt.optMap["USE_GUI"].arg;
	
	// Initializing QT application
	QApplication gemc_gui( argc, argv, (bool) use_gui );

	// Initializing gemc splash class
	// This class will log initialization messages
	// This class will show a splashscreen if use_gui is non zero
	// The screen log verbosity is controlled by LOG_VERBOSITY
	gui_splash gemc_splash(gemcOpt);
	gemc_splash.message(" Initializing GEant4 MonteCarlo");
	
	
	// random seed initialization
	CLHEP::HepRandom::setTheEngine(new CLHEP::MTwistEngine);
	G4int seed;
	
	if(gemcOpt.optMap["RANDOM"].args=="TIME")
	{
		gemc_splash.message(" Initializing CLHEP Random Engine from local time " \
							+ stringify((double) time(NULL)) \
							+ ", cpu clock "        \
							+ stringify((double) clock())    \
							+ " and process id "    \
							+ stringify(getpid()) + ".");
		seed = (G4int) ( (double) time(NULL)- (double) clock()-getpid() );
	}
	else
	{
		seed = atoi(gemcOpt.optMap["RANDOM"].args.c_str());
		gemc_splash.message(" Initializing CLHEP Random Engine from user defined seed.");
	}
	
	CLHEP::HepRandom::setTheSeed(seed);
	gemc_splash.message(" Seed initialized to: " + stringify(seed));
	
	// Construct the default G4 run manager
	gemc_splash.message(" Instantiating Run Manager...");
	G4RunManager *runManager = new G4RunManager;
	
	// Initializing run_condition class
	gemc_splash.message(" Instantiating Run Conditions...");
	runConditions runConds(gemcOpt);
	
	
	// GEMC Detector Map
	gemc_splash.message(" Registering Detectors Factories...");
	// Initializing Detector Factory
	map<string, detectorFactoryInMap> detectorFactoryMap = registerDetectorFactory();
	// Building detector with factories
	map<string, detector> hallMap = buildDetector(detectorFactoryMap, gemcOpt, runConds);
	
	
	// Initialize Materials Map Factory
	gemc_splash.message(" Initializing Material Factories..." );
	map<string, materialFactory> materialFactoriesMap = registerMaterialFactories();
	// Build all materials
	map<string, G4Material*> mats = buildMaterials(materialFactoriesMap, gemcOpt, runConds);
	
	// Initialize Mirrors Map Factory
	gemc_splash.message(" Initializing Mirrors Factories..." );
	map<string, mirrorFactory> mirrorFactoriesMap = registerMirrorFactories();
	// Build all mirrors
	map<string, mirror*> mirs = buildMirrors(mirrorFactoriesMap, gemcOpt, runConds);
	
	
	// Initialize Parameters Map Factory
	gemc_splash.message(" Registering Parameters Factories...");
	map<string, parameterFactoryInMap> parameterFactoriesMap = registerParameterFactories();
	// All Parameters with factories
	map<string, double> gParameters = loadAllParameters(parameterFactoriesMap, gemcOpt, runConds);
	
	
	// Creating the sim_condition map to save to the output
	gemc_splash.message(" Writing simulation parameters in the output...");
	
	// filling gcard option content
	map<string, string> sim_condition = gemcOpt.getOptMap();
	// adding detectors conditions to sim_condition
	mergeMaps(sim_condition, runConds.getDetectorConditionsMap());
	// adding parameters value to sim_condition
	mergeMaps(sim_condition, getParametersMap(gParameters));
	
	// Process Hit Map
	gemc_splash.message(" Building gemc Process Hit Factory...");
	map<string, HitProcess_Factory> hitProcessMap = HitProcess_Map(gemcOpt.optMap["HIT_PROCESS_LIST"].args);

	//addtional hit process
	solid_hitprocess(hitProcessMap);	
	
	///< magnetic Field Map
	gemc_splash.message(" Creating fields Map...");
	map<string, fieldFactoryInMap> fieldFactoryMap = registerFieldFactories();
	map<string, gfield> fieldsMap = loadAllFields(fieldFactoryMap, gemcOpt);
	
	// Build the detector
	gemc_splash.message(" Building Detector Map...");
	MDetectorConstruction* ExpHall = new MDetectorConstruction(gemcOpt);
	ExpHall->hallMap   = &hallMap;
	ExpHall->mirs      = &mirs;
	ExpHall->mats      = &mats;
	ExpHall->fieldsMap = &fieldsMap;
	// this is what calls Construct inside MDetectorConstruction
	runManager->SetUserInitialization(ExpHall);
	
	///< Physics List
	string phys_list = gemcOpt.optMap["PHYSICS"].args  ;
	gemc_splash.message(" Initializing Physics List " + phys_list + "...");
	runManager->SetUserInitialization(new PhysicsList(gemcOpt));

	
	// Setting Max step for all the simulation. This is historically needed to limit
	// the step in magnetic field in vacuum
	double max_step = gemcOpt.optMap["MAX_FIELD_STEP"].arg;
	if(max_step != 0)
		G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetLargestAcceptableStep(max_step);
	
	// Generator
	gemc_splash.message(" Initializing Primary Generator Action...");
	MPrimaryGeneratorAction* gen_action = new MPrimaryGeneratorAction(&gemcOpt);
	runManager->SetUserAction(gen_action);
	
	
	// Event Action
	gemc_splash.message(" Initializing Event Action...");
	MEventAction* event_action = new MEventAction(gemcOpt, gParameters);
	event_action->SetEvtNumber((int) gemcOpt.optMap["EVN"].arg);     ///< Sets event number from OPTION
	runManager->SetUserAction(event_action);
	
	// Stepping Action
	gemc_splash.message(" Initializing Stepping Action...");
	MSteppingAction* SteppingAction = new MSteppingAction(gemcOpt);
	runManager->SetUserAction(SteppingAction);
	
	///< User Interface manager
	gemc_splash.message(" Initializing User Interface...");
	G4UIsession *session = NULL;
	if(use_gui)
		session = new G4UIQt(argc,argv);
	
	///< Vis Manager
	G4VisManager *visManager = NULL;
	if(use_gui)
	{
		visManager = new G4VisExecutive;
		visManager->Initialize();
	}
	
	// Output File: registering output type, output process factory,
	// sensitive detectors into Event Action
	gemc_splash.message(" Initializing Output Action...");
	outputContainer outContainer(gemcOpt);
	map<string, outputFactoryInMap> outputFactoryMap = registerOutputFactories();

	// Initialize G4 kernel
	gemc_splash.message(" Initializing Run Manager...\n");
	// physical volumes, sensitive detectors are built here
	runManager->Initialize();

	// Bank Map, derived from sensitive detector map
	gemc_splash.message(" Creating gemc Banks Map...");
	map<string, gBank> banksMap = read_banks(gemcOpt, runConds.get_systems());
		
	// Getting UI manager, restoring G4Out to cout
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	UImanager->SetCoutDestination(NULL);
	

	// saving simulation condition in the output file
	if(outContainer.outType != "no")
	{
		outputFactory *processOutputFactory  = getOutputFactory(&outputFactoryMap, outContainer.outType);
		processOutputFactory->recordSimConditions(&outContainer, sim_condition);
		// then deleting process output pointer, not needed anymore
		delete processOutputFactory;
	}
	

	event_action->outContainer     = &outContainer;
	event_action->outputFactoryMap = &outputFactoryMap;
	event_action->hitProcessMap    = &hitProcessMap;
	event_action->SeDe_Map         = ExpHall->SeDe_Map;
	event_action->banksMap         = &banksMap;
	event_action->gen_action       = gen_action;
 	
	///< passing output process factory to sensitive detectors
	map<string, sensitiveDetector*>::iterator it;
	for(it = ExpHall->SeDe_Map.begin(); it != ExpHall->SeDe_Map.end(); it++)
		it->second->hitProcessMap = &hitProcessMap;
	
	

	
	gemc_splash.message(" Executing initial directives...\n");
	vector<string> init_commands = init_dmesg(gemcOpt);
	for(unsigned int i=0; i<init_commands.size(); i++)
		UImanager->ApplyCommand(init_commands[i].c_str());
	string exec_macro = "/control/execute " + gemcOpt.optMap["EXEC_MACRO"].args;
	
	clock_t start_events;

	if(use_gui)
	{
		gemc_splash.message("Starting GUI...");
		gemc_gui.processEvents();
		
		gemcMainWidget gemcW(&gemcOpt, runManager, ExpHall->SeDe_Map, &hallMap, mats);
		gemcW.setWindowTitle(GEMC_VERSION);
		gemcW.show();
		
		// splash can finish once gemcW is up
		gemc_splash.splash->finish(&gemcW);
		
		gemc_splash.message(" Executing initial visual directives...\n");
		vector<string> init_vcommands = init_dvmesg(gemcOpt, visManager);
		for(unsigned int i=0; i<init_vcommands.size(); i++)
		{
			gemc_splash.message(" Now executing: " + init_vcommands[i]);
			
			UImanager->ApplyCommand(init_vcommands[i].c_str());
		}
		
		if(exec_macro != "/control/execute no") UImanager->ApplyCommand(exec_macro.c_str());
		if(gemcOpt.optMap["N"].arg>0)
		{
			char command[100];
			sprintf(command, "/run/beamOn %d", (int) gemcOpt.optMap["N"].arg);
			UImanager->ApplyCommand(command);
		}
		
		start_events = clock();
		return gemc_gui.exec();
		// deleting and runManager is now taken care
		// in the gemc_quit slot
		delete visManager;
		if(session != NULL) delete session;
	}
	else
	{
		start_events = clock();
		if(gemcOpt.optMap["N"].arg>0)
		{
			char command[100];
			sprintf(command, "/run/beamOn %d", (int) gemcOpt.optMap["N"].arg);
			UImanager->ApplyCommand(command);
		}
		if(exec_macro != "/control/execute no") UImanager->ApplyCommand(exec_macro.c_str());
	}

	clock_t endTime = clock();
	clock_t clockAllTaken   = endTime - startTime;
	clock_t clockEventTaken = endTime - start_events;

	cout << " > Total gemc time: " <<  clockAllTaken / (double) CLOCKS_PER_SEC << " seconds. "
	     << " Events only time: " << clockEventTaken / (double) CLOCKS_PER_SEC << " seconds. " << endl;

	
	delete runManager;
	return 1;
}








