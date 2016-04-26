void Fun4All_ReadSimExDST(int nEvents=0) {

  std::string dstfile = "simDST.root";
  int runNumber = 429896;
  
  gSystem->Load("libfun4all");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("libbbc.so");
  gSystem->Load("libt0.so");
  gSystem->Load("libmpc.so");
  gSystem->Load("libmpcex_base.so");
//  gSystem->Load("/gpfs/mnt/gpfs02/phenix/mpcex/mrpind/install/lib/libmpcex_base");
  gSystem->Load("libmpcex_interface.so");
//  gSystem->Load("/gpfs/mnt/gpfs02/phenix/mpcex/mrpind/install/lib/libmpcex_modules");
  gSystem->Load("libmpcex_modules.so");
  gSystem->Load("libmpcex_utils.so");
  gSystem->Load("libMyMpcEx.so");
  gSystem->Setenv("ODBCINI","/opt/phenix/etc/odbc.ini.mirror");
  //  gSystem->Load("libjprof.so");
  //  prof *pr = new prof();

  //  getchar();
  
  ///////////////////////////////////////////
  //
  //- Make the Server
  //
  //////////////////////////////////////////
  se = Fun4AllServer::instance();

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER",runNumber);
  rc->set_IntFlag("SIMULATIONFLAG",0x1);
  ///////////////////////////////////////////
  //
  //- Subsystems
  //
  //////////////////////////////////////////

//  rc->set_IntFlag("MPC_RECO_MODE",0x16);
//  MpcReco *mpc = new MpcReco();
//  se->registerSubsystem(mpc);
  
  rc->set_IntFlag("MPCEXFIXEDCALIBS",0x1);
  rc->set_IntFlag("MPCEXCALIBMODE",0x1);
  rc->set_IntFlag("MPCEXCALIBAPPLYSTACK",1);

  mMpcExCreateNodeTree *noder = new mMpcExCreateNodeTree();
  se->registerSubsystem(noder);

  mMpcExDigitizeHits *digitizer = new mMpcExDigitizeHits();
  se->registerSubsystem(digitizer);

  mMpcExSimEventHeader *header = new mMpcExSimEventHeader();
  se->registerSubsystem(header);

  mMpcExLoadCalibrations *calibreader = new mMpcExLoadCalibrations();
  se->registerSubsystem(calibreader);

  SubsysReco *calibdoer = new mMpcExApplyCalibrations();
  se->registerSubsystem(calibdoer);

  mMpcExShower* shower = new mMpcExShower();
  se->registerSubsystem(shower);

//  mMpcExLShower* lshower = new mMpcExLShower();
//  se->registerSubsystem(lshower);

//  mMpcExCreateTree* ct = new mMpcExCreateTree();
//  ct->set_file_name("test.root"); 
//  se->registerSubsystem(ct);
   
    SubsysReco* pi0_reco = new mMpcExPi0Reco();
    se->registerSubsystem(pi0_reco);



  ///////////////////////////////////////////
  //
  //- Analyze the Data.
  //
  //////////////////////////////////////////

  Fun4AllInputManager *simdst = new Fun4AllDstInputManager( "SIM", "DST", "TOP");
//  simdst->AddFile(dstfile);
  se->registerInputManager(simdst);

  string filepath;
//  ifstream myfile("photon_simDst_list.txt");
//  ifstream myfile("simlist.txt");
  ifstream myfile("simlist_pi0.txt");
//  ifstream myfile("simlist_photon.txt");
  int file_counts = 0;
  if(myfile.is_open()){

    gBenchmark->Start("showery");  
 
    while(getline(myfile,filepath)){
      string tfilepath = filepath;
      if(file_counts > 50 ) {
        file_counts++;
	break;
      }
      cout <<"process "<<tfilepath<<endl;
      se->fileopen(simdst->Name(),tfilepath.c_str());
      se->run(0);
      file_counts++;
    }

    se->End();
    gBenchmark->Show("showery");  


    char output[100];
    sprintf(output,"Pi0_mass_sim-%d.root",runNumber);
    Fun4AllHistoManager* hm = se->getHistoManager("Pi0_Mass");
    if(hm) hm->dumpHistos(output);
  }
  else cout<<"file open failed !!!"<<endl;
}
