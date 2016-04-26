// $Id: pisaToDST_IOManager.C,v 1.1 2013/08/06 01:07:18 richi Exp $
/*!
   \file pisaToDST_IOManager.C
   \brief simulation output managers
   \author <a href="mailto:pereira@hep.saclay.cea.fr">Hugo Pereira</a>
   \version $Revision: 1.1 $
   \date $Date: 2013/08/06 01:07:18 $
*/

#include <stdio.h> 
#include <time.h> 

//_______________________________________________________________
void DST_IOManager(const char* dstout, Fun4AllServer* se)
{

  // Control for simulated DST
  Fun4AllDstOutputManager *simDST  = new Fun4AllDstOutputManager("SIMDST", dstout);
 
  recoConsts* rc = recoConsts::instance();

  simDST->AddNode("fkin");
  simDST->AddNode("primary");
  simDST->AddNode("header");

  simDST->AddNode("T0Out");
  simDST->AddNode("VtxOut");
  
  // BBC
  simDST->AddNode("bbcghit");
  simDST->AddNode("BbcOut");
  simDST->AddNode("BbcRaw");

  //MPC
  simDST->AddNode("MpcRaw");

  // DCH
  simDST->AddNode("dcghit");
  simDST->AddNode("DchHitLineTablev1");
  simDST->AddNode("DchHitLineTable");
  simDST->AddNode("DchTrack");
  simDST->AddNode("dDchTracksPerf");
  simDST->AddNode("dDchTracksExtPerf");
  simDST->AddNode("dDchGhitHits");

  // PC
  simDST->AddNode("pc1ghit");
  simDST->AddNode("pc2ghit");
  simDST->AddNode("pc3ghit");
  simDST->AddNode("dPc1GhitClus");
  simDST->AddNode("dPc2GhitClus");
  simDST->AddNode("dPc3GhitClus");
  simDST->AddNode("Pc1Cluster");
  simDST->AddNode("Pc2Cluster");
  simDST->AddNode("Pc3Cluster");
  simDST->AddNode("Pc1Raw");
  simDST->AddNode("Pc2Raw");
  simDST->AddNode("Pc3Raw");

  // RICH
  simDST->AddNode("crkghit");
  simDST->AddNode("CrkHit");
  simDST->AddNode("CrkRing");
  simDST->AddNode("CrkRingBack");

  // TEC
  simDST->AddNode("tecghit");
  simDST->AddNode("dTecGhitRaw");
  simDST->AddNode("TecOut");
  simDST->AddNode("TecHitOut");

  // TOF
  simDST->AddNode("tofghit");
  simDST->AddNode("TofOut");
  simDST->AddNode("dTofGdigi");
  simDST->AddNode("dTofGhitGdigi");
  simDST->AddNode("dTofGdigiRec");
  
  // TOF West
  simDST->AddNode("TofwRaw");
  simDST->AddNode("TofwHit");

  // Aerogel
  simDST->AddNode("AerGeaHits");
  simDST->AddNode("AccCluster");
  simDST->AddNode("AccRaw");
  simDST->AddNode("AccHit");

  // EMCal
  simDST->AddNode("emcghit");
  simDST->AddNode("emcClusterContainer");
  simDST->AddNode("emcTowerContainer");
  simDST->AddNode("emcGeaTrackContainer"); // present only if( EMCSIMULATIONV2 )
  if( 
    rc->FlagExist("EVALUATIONFLAG") && 
    rc->get_IntFlag("EVALUATIONFLAG")==1 )
  { 

    // Evaluation output from EMCal
    simDST->AddNode("dEmcGeaClusterTrack");
    simDST->AddNode("dEmcGeaTrack");
    simDST->AddNode("dEmcGeaTrackCluster");
  }  // check evaluation flag

  // CGL
  simDST->AddNode("CglTrack");
  simDST->AddNode("CglTrackBack");
  simDST->AddNode("PHDchTrackOut");
  simDST->AddNode("PHTrackOut");
  simDST->AddNode("PHTrackOutBack");
  
  // copied from CNT node
  simDST->AddNode("PHCentralTrack");
  simDST->AddNode("PHGlobal");
  simDST->AddNode("PHGlobal_CENTRAL");

  simDST->AddNode("McSingle");

  simDST->AddNode("EventHeader");
  
  se->registerOutputManager(simDST);

  simDST->Print();

  return;

}

//____________________________________________________
void UDST_IOManager(const char* udstout, Fun4AllServer* se)
{

  // Control for simulated microDST (this is not actually used)
  Fun4AllDstOutputManager *usimDST  = new Fun4AllDstOutputManager("USIMDST", udstout);
 
  usimDST->AddNode("TecHitOut");
  usimDST->AddNode("emcTowerContainer");
  usimDST->AddNode("EventHeader");

  se->registerOutputManager(usimDST);

  usimDST->Print();

  return;

}

//____________________________________________________
void CNT_IOManager(const char* ndstout, Fun4AllServer* se)
{

  // Control for simulated CNT
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("SIMCNT", ndstout);
 
  manager->AddNode("PHCentralTrack");
  manager->AddNode("McSingle");
  manager->AddNode("PHGlobal");
  manager->AddNode("PHGlobal_CENTRAL");
  manager->AddNode("EventHeader");
  manager->AddNode("AccCluster"); 
  se->registerOutputManager(manager);

  manager->Print();

  return;

}

//____________________________________________________
void CWG_IOManager(const char* ndstout, Fun4AllServer* se)
{
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("SIMCWG", ndstout);
  manager->AddNode("PhCglList");
  manager->AddNode("PhPhotonList");
  manager->AddNode("HadronPhCglList");
  manager->AddNode("PHGlobal");
  manager->AddNode("PHGlobal_CENTRAL");
  manager->AddNode("EventHeader");
  manager->Print();

  se->registerOutputManager(manager);

  return;
}

//____________________________________________________
void HWG_IOManager(const char* ndstout, Fun4AllServer* se)
{

  // Control for simulated HWG
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("SIMHWG", ndstout);
 
  manager->AddNode("HWGCentralTrack");
  manager->AddNode("McSingle");
  manager->AddNode("PHGlobal");
  manager->AddNode("PHGlobal_CENTRAL");
  manager->AddNode("EventHeader");
  manager->Print();

  se->registerOutputManager(manager);
  return;

}

//____________________________________________________
void EWG_IOManager(const char* ndstout, Fun4AllServer* se)
{

  // Control for simulated HWG
  Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("SIMEWG", ndstout);
 
  manager->AddNode("EWGCentralTrack");
  manager->AddNode("McSingle");
  manager->AddNode("PHGlobal");
  manager->AddNode("PHGlobal_CENTRAL");
  manager->AddNode("RpSumXYObject");
  se->registerOutputManager(manager);

  manager->Print();

  return;

}
