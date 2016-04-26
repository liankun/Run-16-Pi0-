// $Id: pisaToDST.C,v 1.2 2013/12/10 19:41:42 richi Exp $
/*!
   \file pisaToDST.C
   \brief pisa to DST reconstruction chain
   \author <a href="mailto:pereira@hep.saclay.cea.fr">Hugo Pereira</a>
   \version $Revision: 1.2 $
   \date $Date: 2013/12/10 19:41:42 $
*/

void pisaToDST(
  Int_t nEvents = 0, 
  char *filein="PISAEvent.root",
  char *dstout = "simDST.root", 
  int run_number = 429506//last run of Run 13
  
 )
{
 
  // print output
  cout << "pisaToDST - nEvents: " << nEvents << endl;
  cout << "pisaToDST - filein: " << filein << endl;
  if( dstout ) cout << "pisaToDST - dstout: " << dstout << endl;
  cout << "pisaToDST - run_number: " << run_number << endl;
  cout << endl;
   
  // load libraries
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so"); 
  gSystem->Load("librecal.so");
  gSystem->Load("libsimreco.so");
  gSystem->Load("libmpcex_modules.so");
  gSystem->Load("libmuon_util.so" );

  //  getchar();

  ///////////////////////////////////////////
  // recoConsts setup
  //////////////////////////////////////////
    
  recoConsts *rc = recoConsts::instance();

  // 2 means PISA-To-DST
  rc->set_IntFlag("SIMULATIONFLAG",2); 

  // vertex
  rc->set_IntFlag("SIMVERTEXFLAG",2); 
  
  // disable embedding
  rc->set_IntFlag("EMBEDFLAG",0); 
  
  // Reference run number used in 2007 Au+Au 200 GeV 
  rc->set_IntFlag("RUNNUMBER",run_number); 
  
  // Run flag ???
  //  rc->set_IntFlag("RUN8PP200GEV",1);        

  // assume AFS is present as at RCF
  rc->set_IntFlag("AFSABSENT", 0);

  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  ///////////////////////////////////////////
  // Activate the subsystems
  //////////////////////////////////////////
  
  // run header and trigger setting
  se->registerSubsystem( new HeadSimreco() );
  se->registerSubsystem( new TrigSimreco() );

  // BBC simReco
  se->registerSubsystem(new BbcSimreco("BBC"));

  // pisa is used as an input vertex.
  // it overwrites the contents of the BBC out node.
  VtxSimreco* vtx_sim = new VtxSimreco();
  vtx_sim->SmearZ( true ); 
  vtx_sim->UseXY( false );
  vtx_sim->OverwriteBBC( true );  // this is the default 
  vtx_sim->ZVertexSigma(2.); 
  se->registerSubsystem( vtx_sim );
  
  //  TMutExtVtx::get().set_verbosity( MUTOO::SOME );

  // counter
  //MuonCounter* counter = new MuonCounter();
  //  counter->set_event_dump( 100 );
  //  se->registerSubsystem( counter );

  /**
     MPC Reconstruction
     1) Build towers [MPC_RECO_MODE = 0x2] *or 0x3?
     2) Calibrate Simulated Towers 
     3) Build clusters [MPC_RECO_MODE = 0x4
   */
  SubsysReco *mpc_cluster = new MpcReco();

  // 0x4? 0x2?
  rc->set_IntFlag("MPC_RECO_MODE",0x6); //0x1 mpc raw 0x2 mpc tower, 0x4 mpc cluster, 0x7 all
  rc->set_IntFlag("MPC_GAINCORR",0);
  rc->set_IntFlag("MPC_EVENT_DISPLAY",0);
  rc->set_IntFlag("MPC_VERBOSITY",0);

  se->registerSubsystem(mpc_cluster);

  se->registerSubsystem(new VtxReco("VTX"));
  se->registerSubsystem(new T0Reco()); 

  //  This is the class which makes the GlobalEvent data on the nanoDST output
  se->registerSubsystem(new GlobalReco());
  se->registerSubsystem(new NCCSimreco("NCC"));

  mMpcExCreateNodeTree *nodes = new mMpcExCreateNodeTree();
  se->registerSubsystem(nodes);
  mMpcExUnpackPISA *unpacker = new mMpcExUnpackPISA();
  se->registerSubsystem(unpacker);

  ///////////////////////////////////////////
  // InputManager 
  ///////////////////////////////////////////
  Fun4AllInputManager *inMan = new Fun4AllPisaInputManager("PisaIn","TOP");
  se->registerInputManager(inMan);

  ///////////////////////////////////////////
  // open input file
  inMan->AddFile(filein); // load the filelist into the Input Manager

  ///////////////////////////////////////////
  // OutputManagers Set up functions  
  ///////////////////////////////////////////
  
  Fun4AllDstOutputManager *simDST  = new Fun4AllDstOutputManager("SIMDST", dstout);
  se->registerOutputManager(simDST);

  //add fkin and primary


  simDST->AddNode("Sync");
  simDST->AddNode("RunHeader");
  simDST->AddNode("EventHeader");
   //add fkin and primary

  simDST->AddNode("fkin");
  simDST->AddNode("primary");
  
  simDST->AddNode("T0Out");
  simDST->AddNode("VtxOut");
  simDST->AddNode("BbcOut");
  simDST->AddNode("BbcRaw");
  simDST->AddNode("BbcGeo");
  simDST->AddNode("PHGlobal");
  simDST->AddNode("mpcClusterContainer");
  simDST->AddNode("mpcTowerContainer");
  simDST->AddNode("MpcRaw");
  simDST->AddNode("TrigLvl1");
  simDST->AddNode("TMpcExGeaHitContainer"); //These are the Geant Energy Hits

  gBenchmark->Start("eventLoop");
  // process input events
  se->run(nEvents);
  se->End();
  gBenchmark->Show("eventLoop");

  // If you do not see this message, the job failed
  cout << "Completed reconstruction." << endl;
}
