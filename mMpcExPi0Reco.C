#include "mMpcExPi0Reco.h"
#include "MpcExConstants.h"
#include "TMpcExHitContainer.h"
#include "TMpcExHit.h"
#include "TMpcExLShower.h"
#include "TMpcExLShowerContainer.h"
#include "TMpcExShowerContainer.h"
#include <MpcExMapper.h>
#include <PHIODataNode.h> 
#include <getClass.h>
#include <BbcOut.h>
#include <PHGlobal.h>
#include <Bbc.hh>
#include <TriggerHelper.h>
#include <Exogram.h>

#include <MpcMap.h>
#include <mpcClusterContent.h>
#include <mpcClusterContainer.h>
#include <mpcClusterContentV1.h>
#include <mpcTowerContainer.h>
#include <mpcTowerContent.h>
#include <mpcTowerContentV1.h>

#include <primary.h>
#include <primaryWrapper.h>
#include <fkinWrapper.h>

#include <Fun4AllReturnCodes.h>
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <recoConsts.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>

using namespace std;
using namespace findNode;

mMpcExPi0Reco::mMpcExPi0Reco(const char* name) :
  SubsysReco(name){
  _vertex = -9999.9;
  _mpcex_hit_container = NULL;
  _mpcex_shower_container = NULL;
  _mpcex_lshower_container = NULL;
  _mpc_cluster_container = NULL;
  _mpc_tower_container = NULL;
  _mpc_map = NULL;
  for(int iarm = 0;iarm < 2;iarm++){
    grammyh[iarm] = NULL;
    grammyl[iarm] = NULL;
    hmpc_gridxy[iarm] = NULL;
  }

  }

int mMpcExPi0Reco::End(PHCompositeNode* topNode){
  return EVENT_OK;
}

mMpcExPi0Reco::~mMpcExPi0Reco(){

}

int mMpcExPi0Reco::Init(PHCompositeNode* topNode){
  pi0_mass_sim_init();

  return EVENT_OK;
}

int mMpcExPi0Reco::InitRun(PHCompositeNode* topNode){
  set_interface_ptrs(topNode);
  return EVENT_OK;
}

void mMpcExPi0Reco::set_interface_ptrs(PHCompositeNode* topNode){


   _mpcex_hit_container = getClass<TMpcExHitContainer>(topNode,"TMpcExHitContainer");
  if(!_mpcex_hit_container){
    cout << PHWHERE <<":: No TMpcExHitContainer!!!"<<endl;
    exit(1);
  }


  _mpcex_shower_container = getClass<TMpcExShowerContainer>(topNode,"TMpcExShowerContainer");
  
  if(!_mpcex_shower_container){
    cout << PHWHERE <<":: No TMpcExShowerContainer !!!"<<endl;
    exit(1);
  }

  _mpcex_lshower_container = getClass<TMpcExLShowerContainer>(topNode,"TMpcExLShowerContainer");
  if(!_mpcex_lshower_container){
    cout <<PHWHERE <<":: No TMpcExLShowerContainer !!!"<<endl;    
//    exit(1);
  }

  _mpc_cluster_container = getClass<mpcClusterContainer>(topNode,"mpcClusterContainer");
  if(!_mpc_cluster_container){
    cout <<PHWHERE <<":: No mpcClusterContainer!!!"<<endl;
    exit(1);
  }
  _mpc_tower_container = getClass<mpcTowerContainer>(topNode,"mpcTowerContainer");
  if(!_mpc_tower_container){
    cout << PHWHERE <<":: No mpcTowerContainer!!!"<<endl;
    exit(1);
  }
  _mpc_map = MpcMap::instance();
  if(!_mpc_map){
    cout <<"No MpcMap!!!"<<endl;
    exit(1);
  }
}


int mMpcExPi0Reco::process_event(PHCompositeNode* topNode){
  PHGlobal* phglobal = getClass<PHGlobal>(topNode,"PHGlobal");
  BbcOut* bbcout = getClass<BbcOut>(topNode,"BbcOut");
  if(!bbcout && !phglobal){
    cout <<"No BbcOut or PHGlobal !!!"<<endl;
    exit(1);
  }
  _vertex = (phglobal==0) ? phglobal->getBbcZVertex() : bbcout->get_VertexPoint();

//  cout <<"Pi Reco "<<endl;
//  pi0_mass_sim(topNode);
//  pi0_mass_sim_PP(topNode);
//  pi0_mass_data(topNode);
//  event_display(topNode);
  pi0_mass_simV2_study(topNode);
  return EVENT_OK;
}

int mMpcExPi0Reco::pi0_mass_sim_init(){
  Fun4AllServer* se = Fun4AllServer::instance();
  Fun4AllHistoManager* hm = se->getHistoManager("Pi0_Mass");
  if(!hm){
    hm = new Fun4AllHistoManager("Pi0_Mass");
    se->registerHistoManager(hm);
  }



  hPi0_mass[0] = new TH1D("hPi0_mass_arm0","Pi0 Mass Arm 0",500,0,0.5);
  hPi0_mass[0]->GetXaxis()->SetTitle("Mass/GeV");
  hm->registerHisto(hPi0_mass[0]);

  hPi0_mass[1] = new TH1D("hPi0_mass_arm1","Pi0 Mass Arm 1",500,0,0.5);
  hPi0_mass[1]->GetYaxis()->SetTitle("Mass/GeV");
  hm->registerHisto(hPi0_mass[1]);
/*
  hAsy_tote_angle = new TH3D("hAsy_tote_angle","Asymmetry Energy Agnle",100,0,1,200,0,100,100,0,0.1);
  hAsy_tote_angle->GetXaxis()->SetTitle("Asymmetry");
  hAsy_tote_angle->GetYaxis()->SetTitle("Energy/GeV");
  hAsy_tote_angle->GetZaxis()->SetTitle("dhr");
  hm->registerHisto(hAsy_tote_angle);

  hmass_dr[0] = new TH2D("hmass_dr_arm0","Mass vs dr Arm 0",500,0,0.5,500,0,0.05);
  hmass_dr[0]->GetXaxis()->SetTitle("Mass/GeV");
  hmass_dr[0]->GetYaxis()->SetTitle("hsdr");
  hm->registerHisto(hmass_dr[0]);

  hmass_dr[1] = new TH2D("hmass_dr_arm1","Mass vs dr Arm 1",500,0,0.5,500,0,0.05);
  hmass_dr[1]->GetXaxis()->SetTitle("Mass/GeV");
  hmass_dr[1]->GetYaxis()->SetTitle("hsdr");
  hm->registerHisto(hmass_dr[1]);

  hmass_angle[0] = new TH2D("hmass_angle_arm0","Mass vs Open Angle Arm 0",500,0,0.5,100,0,0.05);
  hmass_angle[0]->GetXaxis()->SetTitle("Mass/GeV");
  hmass_angle[0]->GetYaxis()->SetTitle("Open Angle in hsdr");
  hm->registerHisto(hmass_angle[0]);

  hmass_angle[1] = new TH2D("hmass_angle_arm1","Mass vs Open Angle Arm 1",500,0,0.5,100,0,0.05);
  hmass_angle[1]->GetXaxis()->SetTitle("Mass/GeV");
  hmass_angle[1]->GetYaxis()->SetTitle("Open Angle in hsdr");
  hm->registerHisto(hmass_angle[1]);

  he_reco_real0[0] = new TH2D("he_reco_real0_arm0","first photon Reco E vs True E arm 0",100,0,30,100,0,30);
  he_reco_real0[0]->GetXaxis()->SetTitle("Reco E/GeV");
  he_reco_real0[0]->GetYaxis()->SetTitle("Real E/GeV");
  hm->registerHisto(he_reco_real0[0]);

  he_reco_real0[1] = new TH2D("he_reco_real0_arm1","first photon Reco E vs True E arm 1",100,0,30,100,0,30);
  he_reco_real0[1]->GetXaxis()->SetTitle("Reco E/GeV");
  he_reco_real0[1]->GetYaxis()->SetTitle("Real E/GeV");
  hm->registerHisto(he_reco_real0[1]);

  he_reco_real1[0] = new TH2D("he_reco_real1_arm0","second photon Reco E vs True E arm 0",100,0,30,100,0,30);
  he_reco_real0[0]->GetXaxis()->SetTitle("Reco E/GeV");
  he_reco_real0[0]->GetYaxis()->SetTitle("Real E/GeV");
  hm->registerHisto(he_reco_real1[0]);

  he_reco_real1[1] = new TH2D("he_reco_real1_arm1","second photon Reco E vs True E arm 1",100,0,30,100,0,30);
  he_reco_real1[1]->GetXaxis()->SetTitle("Reco E/GeV");
  he_reco_real1[1]->GetYaxis()->SetTitle("Real E/GeV");
  hm->registerHisto(he_reco_real1[1]);

  hnomatch_e_ratio[0] = new TH2D("hnomatch_e_ratio_arm0","Nomatch ratio Arm 0",100,0,5,100,0,5);
  hnomatch_e_ratio[0]->GetXaxis()->SetTitle("NoMatch ratio");
  hnomatch_e_ratio[0]->GetYaxis()->SetTitle("Match ratio");
  hm->registerHisto(hnomatch_e_ratio[0]);

  hnomatch_e_ratio[1] = new TH2D("hnomatch_e_ratio_arm1","Nomatch ratio Arm 1",100,0,5,100,0,5);
  hnomatch_e_ratio[1]->GetXaxis()->SetTitle("NoMatch ratio");
  hnomatch_e_ratio[1]->GetYaxis()->SetTitle("Match ratio");
  hm->registerHisto(hnomatch_e_ratio[1]);

  hnomatch_rms[0] = new TH2D("hnomatch_rms_arm0","NoMatch rms Arm 0",200,0,0.02,200,0,0.02);
  hnomatch_rms[0]->GetXaxis()->SetTitle("NoMatch RMS");
  hnomatch_rms[0]->GetYaxis()->SetTitle("Match RMS");
  hm->registerHisto(hnomatch_rms[0]);

  hnomatch_rms[1] = new TH2D("hnomatch_rms_arm1","NoMatch rms Arm 1",200,0,0.02,200,0,0.02);
  hnomatch_rms[1]->GetXaxis()->SetTitle("NoMatch RMS");
  hnomatch_rms[1]->GetYaxis()->SetTitle("Match RMS");
  hm->registerHisto(hnomatch_rms[1]);

  hnomatch_cluster[0] = new TH2D("hnomatch_cluster_arm0","NoMatch cluster Arm 0",200,0,0.2,200,0,0.2);
  hnomatch_cluster[0]->GetXaxis()->SetTitle("NoMatch");
  hnomatch_cluster[0]->GetYaxis()->SetTitle("Match");
  hm->registerHisto(hnomatch_cluster[0]);

  hnomatch_cluster[1] = new TH2D("hnomatch_cluster_arm1","NoMatch cluster Arm 1",200,0,0.2,200,0,0.2);
  hnomatch_cluster[1]->GetXaxis()->SetTitle("NoMatch");
  hnomatch_cluster[1]->GetYaxis()->SetTitle("Match");
  hm->registerHisto(hnomatch_cluster[1]);

  hnomatch_fired_layers = new TH2D("hnomatch_fired_layers","Fired layers",9,-0.5,8.5,9,-0.5,8.5);
  hnomatch_fired_layers->GetXaxis()->SetTitle("nomatch");
  hnomatch_fired_layers->GetYaxis()->SetTitle("match");
  hm->registerHisto(hnomatch_fired_layers);

  hnomatch_Nhits = new TH2D("hnomatch_Nhits","Number of hits",100,0,500,100,0,500);
  hnomatch_Nhits->GetXaxis()->SetTitle("nomatch");
  hnomatch_Nhits->GetYaxis()->SetTitle("match");
  hm->registerHisto(hnomatch_Nhits);

  htote_reco_real = new TH2D("htote_reco_real","Pi0 Energy Reco vs Real",200,0,100,200,0,100);
  htote_reco_real->GetXaxis()->SetTitle("Real E/GeV");
  htote_reco_real->GetYaxis()->SetTitle("Reconstruction E/GeV");
  hm->registerHisto(htote_reco_real);

  hangle_reco_real = new TH2D("hangle_reco_real","Open Angle Reco vs Real",100,0,0.5,100,0,0.5);
  hangle_reco_real->GetXaxis()->SetTitle("Real Angle");
  hangle_reco_real->GetYaxis()->SetTitle("Reco Angle");
  hm->registerHisto(hangle_reco_real);

  htote_clus = new TH2D("htote_clus","Pi0 Energy Reco vs Real (merge shower)",200,0,100,200,0,100);
  htote_clus->GetXaxis()->SetTitle("Real E/GeV");
  htote_clus->GetYaxis()->SetTitle("Reconstruction E/GeV");
  hm->registerHisto(htote_clus);

  hangle_clus = new TH2D("hangle_clus","Open Angle Reco vs Real (merge shower)",100,0,0.5,100,0,0.5);
  hangle_clus->GetXaxis()->SetTitle("Real Angle");
  hangle_clus->GetYaxis()->SetTitle("Reco Angle");
  hm->registerHisto(hangle_clus);
 */
  return EVENT_OK;
}

int mMpcExPi0Reco::pi0_mass_sim(PHCompositeNode* topNode){
 
  fkinWrapper* fkin = findNode::getClass<fkinWrapper>(topNode, "fkin");
  primaryWrapper* primary = findNode::getClass<primaryWrapper>(topNode, "primary");
  
  
  if(!fkin || !primary){
    cout <<PHWHERE<<"No fkinWrapper !!!"<<endl;
    return EVENT_OK;
  }
  size_t fkinrows = fkin->RowCount();
  size_t primrows = fkin->RowCount();
  //it looks phi and theta is in degree
  double tot_e0 = 0;
  double theta0 = 0;
  double phi0 = 0;
  double tot_e1 = 0;
  double theta1 = 0;
  double phi1 = 0;
  bool isfirst = true;
  bool isPi0 = false;

  
  //no decay photon is recorded in primary
//  for(size_t iprim = 0;iprim < primrows;iprim++){
//    int idpart = primary->get_idpart(iprim);
//    cout <<idpart<<endl;
//    if(iprim > 3) break;
//  }
  
  for(size_t ifkin = 0;ifkin < fkinrows;ifkin++){
    int idpart = fkin->get_idpart(ifkin);
    int idparent = fkin->get_idparent(ifkin);
//    cout <<idpart<<" "<<idparent<<endl;
    //get the photon
    if(idpart != 1) continue;
    //pi0 
    if(idparent == 7){
//      cout <<"find pi0 !"<<endl;
      if(isfirst){
        tot_e0 = fkin->get_ptot(ifkin);
	theta0 = fkin->get_pthet(ifkin)*3.14159/180.;
	phi0 = fkin->get_pphi(ifkin)*3.14159/180.;
//	cout <<"first photon: "<<tot_e0<<" "<<theta0<<" "<<phi0<<endl;
	isfirst = false;
      }
      else{
        tot_e1 = fkin->get_ptot(ifkin);
	theta1 = fkin->get_pthet(ifkin)*3.14159/180.;
	phi1 = fkin->get_pphi(ifkin)*3.14159/180.;
	isPi0 = true;
//        cout <<"second photon: "<<tot_e1<<" "<<theta1<<" "<<phi1<<endl;
	break;
      }
    }
  }

  double asy = (tot_e0 - tot_e1)/(tot_e0 + tot_e1);
  double tot_e = tot_e0 + tot_e1;
  double rhsx0 = -sin(theta0)*cos(phi0)/cos(theta0);
  double rhsy0 = sin(theta0)*sin(phi0)/cos(theta0);

  double rhsx1 = -sin(theta1)*cos(phi1)/cos(theta1);
  double rhsy1 = sin(theta1)*sin(phi1)/cos(theta1);

//  cout <<"photon 0: "<<tot_e0<<" "<<rhsx0<<" "<<rhsy0<<" "<<cos(theta0)<<endl;
//  cout <<"photon 1: "<<tot_e1<<" "<<rhsx1<<" "<<rhsy1<<" "<<cos(theta1)<<endl;
  
  double hdr = (rhsx0-rhsx1)*(rhsx0-rhsx1)+(rhsy0-rhsy1)*(rhsy0-rhsy1);
  hdr = sqrt(hdr);
  
  hAsy_tote_angle->Fill(asy,tot_e,hdr);

/*************************
  //test algebra : the algebra is right
  if(!isPi0) return EVENT_OK;
  double hsx0 = sin(theta0)*cos(phi0)/cos(theta0);
  double hsy0 = sin(theta0)*sin(phi0)/cos(theta0);

  double hsx1 = sin(theta1)*cos(phi1)/cos(theta1);
  double hsy1 = sin(theta1)*sin(phi1)/cos(theta1);
  double norm0 = sqrt(hsx0*hsx0+hsy0*hsy0+1*1);
  double norm1 = sqrt(hsx1*hsx1+hsy1*hsy1+1*1);
  double prdct = hsx0*hsx1+hsy0*hsy1+1*1;
  double cs = prdct/(norm0*norm1);
  double mass = sqrt(2*tot_e0*tot_e1*(1-cs));
//  cout <<"mass: "<<mass<<" cos "<<cs<<endl;
  hPi0_mass[0]->Fill(mass);
*****************************/

  if(tot_e0+tot_e1 > 25 ) return EVENT_OK;



  int Nshowers = _mpcex_shower_container->size();
  if(Nshowers < 2 || Nshowers > 100) return EVENT_OK;
//  cout <<"there are "<<Nshowers<<endl;
  int match0 = -1;
  int match1 = -1;
  double smt_dr0 = 9999;
  for(int i = 0;i < Nshowers;i++){
    TMpcExShower* shower = _mpcex_shower_container->getShower(i);
    double hsx = shower->get_hsx();
    double hsy = shower->get_hsy();
    if(shower->get_CalibEInRange() == 0)continue;
    double dist = (hsx-rhsx0)*(hsx-rhsx0)+(hsy-rhsy0)*(hsy-rhsy0); 
    dist = sqrt(dist);
    if(smt_dr0 > dist){
      smt_dr0 = dist;
      match0 = i;
    }
  }
//  cout <<smt_dr0<<endl;  
  double smt_dr1 = 9999;
  for(int i = 0;i < Nshowers;i++){
    if(i == match0) continue;
    TMpcExShower* shower = _mpcex_shower_container->getShower(i);
    double hsx = shower->get_hsx();
    double hsy = shower->get_hsy();
    double dist = (hsx-rhsx1)*(hsx-rhsx1)+(hsy-rhsy1)*(hsy-rhsy1); 
    if(shower->get_CalibEInRange() == 0)continue;
    dist = sqrt(dist);
    if(smt_dr1 > dist){
      smt_dr1 = dist;
      match1 = i;
    }
  }
//  cout<<smt_dr1<<endl;
  
  bool good_evt = false;
  double match_ratio0 = 0;
  double rms0 = 0;
  double clus_dist0 = 0;
  double fired_layers0 = 0;
  double Nhits0 = 0;
  if(match1!= -1 && match0!=-1){
    TMpcExShower* match_shower0 = _mpcex_shower_container->getShower(match0);
    double match_e0 = match_shower0->get_roughTotE();
    double match_hsx0 = match_shower0->get_hsx();
    double match_hsy0 = match_shower0->get_hsy();
    int match_arm0 = match_shower0->get_arm();
    int match_calib0 = match_shower0->get_CalibEInRange();
    int match_clus0 = match_shower0->get_ClosestMPCClusterIndex();

    match_ratio0 = match_shower0->get_esum()/match_shower0->get_mpcECorr();
    
    double rms_x0 = match_shower0->get_rms_hsx();
    double rms_y0 = match_shower0->get_rms_hsy();
    rms0 = sqrt(rms_x0*rms_x0+rms_y0*rms_y0);
    clus_dist0 = match_shower0->get_ClosestMPCClusterDistance();
    fired_layers0 = match_shower0->get_nlayers();
    Nhits0 = match_shower0->sizeHits();

    TMpcExShower* match_shower1 = _mpcex_shower_container->getShower(match1);
    double match_e1 = match_shower1->get_roughTotE();
    double match_hsx1 = match_shower1->get_hsx();
    double match_hsy1 = match_shower1->get_hsy();
    int match_arm1 = match_shower1->get_arm();
    int match_calib1 = match_shower1->get_CalibEInRange();
    int match_clus1 = match_shower1->get_ClosestMPCClusterIndex();

    
    double match_hsdr = (match_hsx0-match_hsx1)*(match_hsx0-match_hsx1)+(match_hsy0-match_hsy1)*(match_hsy0-match_hsy1);
    match_hsdr = sqrt(match_hsdr);

    double match_norm0 = sqrt(match_hsx0*match_hsx0+match_hsy0*match_hsy0+1*1);
    double match_norm1 = sqrt(match_hsx1*match_hsx1+match_hsy1*match_hsy1+1*1);
    double match_prdct = match_hsx0*match_hsx1+match_hsy0*match_hsy1+1*1;
    double match_cs = match_prdct/(match_norm0*match_norm1);
    double match_mass = sqrt(2*match_e0*match_e1*(1-match_cs));    

    //use correct open angle test the energy calibration


    //smt_dr to cuts not match photons
    if(match_arm0 ==  match_arm1 && match_calib0 && match_calib1 && (smt_dr0+smt_dr1)/2. < 0.005){
      good_evt = true;
      hmass_dr[match_arm0]->Fill(match_mass,(smt_dr0+smt_dr1)/2.);
      hmass_angle[match_arm0]->Fill(match_mass,hdr);
      he_reco_real0[match_arm0]->Fill(match_e0,tot_e0);
      he_reco_real1[match_arm0]->Fill(match_e1,tot_e1);
      htote_reco_real->Fill(tot_e1+tot_e0,match_e0+match_e1);
      hangle_reco_real->Fill(hdr,match_hsdr);
//      cout <<(smt_dr0+smt_dr1)/2.<<endl;
      //merge the shower if they have the same cluster
     
      double merge_hsx0 = match_hsx0*match_e0;
      double merge_hsy0 = match_hsy0*match_e0;
      double merge_mpcex_e0 = match_shower0->get_esum();

      double merge_hsx1 = match_hsx1*match_e1;
      double merge_hsy1 = match_hsy1*match_e1;
      double merge_mpcex_e1 = match_shower1->get_esum();


      for(int i2 = 0;i2 < Nshowers;i2++){
        TMpcExShower* myshower = _mpcex_shower_container->getShower(i2);
	int clus_index = myshower->get_ClosestMPCClusterIndex();
	if(clus_index == match_clus0){
          merge_hsx0 += myshower->get_hsx()*myshower->get_esum();
	  merge_hsy0 += myshower->get_hsy()*myshower->get_esum();
	  merge_mpcex_e0 += myshower->get_esum();
	}
	if(clus_index == match_clus1){
          merge_hsx1 += myshower->get_hsx()*myshower->get_esum();
	  merge_hsy1 += myshower->get_hsy()*myshower->get_esum();
	  merge_mpcex_e1 += myshower->get_esum();
	}
      }//i2
      if(match_clus0 != match_clus1){
        merge_hsx0 = merge_hsx0/merge_mpcex_e0;
	merge_hsy0 = merge_hsy0/merge_mpcex_e0;
	merge_hsx1 = merge_hsx1/merge_mpcex_e1;
	merge_hsy1 = merge_hsy1/merge_mpcex_e1;
	mpcClusterContent* clus0 = _mpc_cluster_container->getCluster(match_clus0);
        mpcClusterContent* clus1 = _mpc_cluster_container->getCluster(match_clus1);

        double merge_dr = (merge_hsx0-merge_hsx1)*(merge_hsx0-merge_hsx1)+(merge_hsy0-merge_hsy1)*(merge_hsy0-merge_hsy1);
	merge_dr = sqrt(merge_dr);
        htote_clus->Fill(tot_e0+tot_e1,clus0->e()+clus1->e()+merge_mpcex_e0+merge_mpcex_e1);
	hangle_clus->Fill(hdr,merge_dr);
      }
    
    }
//    else cout <<"Not Match !!!"<<endl;
  }

  if(!good_evt) return EVENT_OK;
  
  for(int i = 0;i < Nshowers;i++){
    TMpcExShower* shower0 = _mpcex_shower_container->getShower(i);
    
    if(shower0->get_CalibEInRange() == 0) continue;
    //esum is the corrected MpcEx energy
    double e0 = shower0->get_roughTotE();
    
//    int clus_index0 = shower0->get_ClosestMPCClusterIndex();
//    cout <<clus_index0<<endl;
//    mpcClusterContent* clus0 = _mpc_cluster_container->getCluster(clus_index0);
//    double clus0_e = clus0->e();
//    e0 = shower0->get_esum()+clus0_e;
    if(e0 < 0.5 || e0 > 25) continue;
    double hsx0 = shower0->get_hsx();
    double hsy0 = shower0->get_hsy();
    int arm0 = shower0->get_arm();
    double mpcE3x30 = shower0->get_mpcE3x3();
    double mpcE5x50 = shower0->get_mpcE5x5();
    
    double ratio = shower0->get_esum()/shower0->get_mpcECorr();
    double rms_x = shower0->get_rms_hsx();
    double rms_y = shower0->get_rms_hsy();
    double rms = sqrt(rms_x*rms_x+rms_y*rms_y);
    double clus_dist = shower0->get_ClosestMPCClusterDistance();
    
    if(i != match0 && i != match1){
      hnomatch_e_ratio[arm0]->Fill(ratio,match_ratio0);
      hnomatch_rms[arm0]->Fill(rms,rms0);
      hnomatch_cluster[arm0]->Fill(clus_dist,clus_dist0);
      hnomatch_fired_layers->Fill(shower0->get_nlayers(),fired_layers0);
      hnomatch_Nhits->Fill(shower0->sizeHits(),Nhits0);
    }
    if(shower0->get_nlayers() < 5) continue;
    if(shower0->sizeHits() < 15) continue;
    if(ratio < 0.015 || ratio > 0.9 || rms < 0.0015 || rms > 0.004|| clus_dist > 0.03|| (mpcE3x30/mpcE5x50 < 0.9)) continue;

//    cout <<arm0<<" "<<e0<<" "<<hsx0<<" "<<hsy0<<endl;
    
    for(int j = i+1;j < Nshowers;j++){
     
      TMpcExShower* shower1 = _mpcex_shower_container->getShower(j);
      if(shower1->get_CalibEInRange() == 0) continue;
      double e1 = shower1->get_roughTotE();
//      int clus_index1 = shower1->get_ClosestMPCClusterIndex();
//      cout << clus_index1<<endl;
//      mpcClusterContent* clus1 = _mpc_cluster_container->getCluster(clus_index1);
//      double clus1_e = clus1->e();
//      e1 = shower1->get_esum()+clus1_e;
     
      if(e1 < 0.5 || e1 > 25) continue;
      int arm1 = shower1->get_arm();
      if(arm0 != arm1) continue;
      if(shower1->get_nlayers() < 5) continue;
      if(shower1->sizeHits() < 15) continue;

      double hsx1 = shower1->get_hsx();
      double hsy1 = shower1->get_hsy();
      double mpcE3x31 = shower1->get_mpcE3x3();
      double mpcE5x51 = shower1->get_mpcE5x5();
//      double hsdr = (hsx1-hsx0)*(hsx1-hsx0)+(hsy1-hsy0)*(hsy1-hsy0);
//      hsdr = sqrt(hsdr);
//      if(hsdr < 0.01) continue;
      ratio = shower1->get_esum()/shower1->get_mpcECorr();
      rms_x = shower1->get_rms_hsx();
      rms_y = shower1->get_rms_hsy();
      rms = sqrt(rms_x*rms_x+rms_y*rms_y);
      clus_dist = shower1->get_ClosestMPCClusterDistance();
      if(ratio < 0.015 || ratio > 0.9 || rms < 0.0015 || rms > 0.004|| clus_dist > 0.03 || (mpcE5x51/mpcE3x31 < 0.9)) continue;

      double norm0 = sqrt(hsx0*hsx0+hsy0*hsy0+1*1);
      double norm1 = sqrt(hsx1*hsx1+hsy1*hsy1+1*1);
      double prdct = hsx0*hsx1+hsy0*hsy1+1*1;
      double cs = prdct/(norm0*norm1);
      double mass = sqrt(2*e0*e1*(1-cs));
//      cout <<"mass: "<<mass<<endl;
      hPi0_mass[arm0]->Fill(mass);
    }
  }
  
  return EVENT_OK;  
}

int mMpcExPi0Reco::pi0_mass_sim_PP(PHCompositeNode* topNode){
  int Nshowers = _mpcex_shower_container->size();
//  if(Nshowers < 2 || Nshowers > 8) return EVENT_OK;

  for(int i = 0;i < Nshowers;i++){
    TMpcExShower* shower0 = _mpcex_shower_container->getShower(i);
    if(shower0->get_CalibEInRange() == 0) continue;
    if(shower0->get_nlayers() < 3) continue;
    if(shower0->sizeHits() < 3) continue;
    //esum is the corrected MpcEx energy
    double e0 = shower0->get_roughTotE();
    if(e0 < 0.5 || e0 > 25) continue;
    double hsx0 = shower0->get_hsx();
    double hsy0 = shower0->get_hsy();
    int arm0 = shower0->get_arm();
    double ratio = shower0->get_esum()/shower0->get_mpcECorr();
    double rms_x = shower0->get_rms_hsx();
    double rms_y = shower0->get_rms_hsy();
    double rms = sqrt(rms_x*rms_x+rms_y*rms_y);
    double clus_dist = shower0->get_ClosestMPCClusterDistance();

    if(ratio < 0.015 || ratio > 0.9 || rms < 0.0015 || rms > 0.004|| clus_dist > 0.015) continue;


//    cout <<arm0<<" "<<e0<<" "<<hsx0<<" "<<hsy0<<endl;
    for(int j = i+1;j < Nshowers;j++){
      TMpcExShower* shower1 = _mpcex_shower_container->getShower(j);
      if(shower1->get_CalibEInRange() == 0) continue;
      double e1 = shower1->get_roughTotE();
      if(e1 < 0.5 || e1 > 25 || e0+e1 > 25) continue;
      int arm1 = shower1->get_arm();
      if(arm0 != arm1) continue;
      if(shower1->get_nlayers() < 3) continue;
      if(shower1->sizeHits() < 3) continue;

      double hsx1 = shower1->get_hsx();
      double hsy1 = shower1->get_hsy();
//      double hsdr = (hsx1-hsx0)*(hsx1-hsx0)+(hsy1-hsy0)*(hsy1-hsy0);
//      hsdr = sqrt(hsdr);
//      if(hsdr < 0.01) continue;
      ratio = shower1->get_esum()/shower1->get_mpcECorr();
      rms_x = shower1->get_rms_hsx();
      rms_y = shower1->get_rms_hsy();
      rms = sqrt(rms_x*rms_x+rms_y*rms_y);
      clus_dist = shower1->get_ClosestMPCClusterDistance();
      if(ratio < 0.015 || ratio > 0.9 || rms < 0.0015 || rms > 0.004|| clus_dist > 0.015) continue;

      double norm0 = sqrt(hsx0*hsx0+hsy0*hsy0+1*1);
      double norm1 = sqrt(hsx1*hsx1+hsy1*hsy1+1*1);
      double prdct = hsx0*hsx1+hsy0*hsy1+1*1;
      double cs = prdct/(norm0*norm1);
      double mass = sqrt(2*e0*e1*(1-cs));
//      cout <<"mass: "<<mass<<endl;
      hPi0_mass[arm0]->Fill(mass);
    }
  }

  return 0;
}

int mMpcExPi0Reco::pi0_mass_data(PHCompositeNode* topNode){
  int Nshowers = _mpcex_shower_container->size();
  if(Nshowers < 2) return EVENT_OK;
  cout << "Number of showers: "<<Nshowers<<endl;
 
 
  for(int i = 0;i < Nshowers;i++){
    TMpcExShower* shower0 = _mpcex_shower_container->getShower(i);
    if(shower0->get_CalibEInRange() == 0) continue;
    if(shower0->get_nlayers() < 3) continue;
    if(shower0->sizeHits() < 3) continue;
    //esum is the corrected MpcEx energy
    double e0 = shower0->get_roughTotE();
    if(e0 < 0.5 || e0 > 25) continue;
    double hsx0 = shower0->get_hsx();
    double hsy0 = shower0->get_hsy();
    int arm0 = shower0->get_arm();
    double ratio = shower0->get_esum()/shower0->get_mpcECorr();
    double rms_x = shower0->get_rms_hsx();
    double rms_y = shower0->get_rms_hsy();
    double rms = sqrt(rms_x*rms_x+rms_y*rms_y);
    double clus_dist = shower0->get_ClosestMPCClusterDistance();

    if(ratio < 0.015 || ratio > 0.9 || rms < 0.0015 || rms > 0.004|| clus_dist > 0.003) continue;


//    cout <<arm0<<" "<<e0<<" "<<hsx0<<" "<<hsy0<<endl;
    for(int j = i+1;j < Nshowers;j++){
      TMpcExShower* shower1 = _mpcex_shower_container->getShower(j);
      if(shower1->get_CalibEInRange() == 0) continue;
      double e1 = shower1->get_roughTotE();
      if(e1 < 0.5 || e1 > 25 || e0+e1 > 25) continue;
      int arm1 = shower1->get_arm();
      if(arm0 != arm1) continue;
      if(shower1->get_nlayers() < 4) continue;
      if(shower1->sizeHits() < 10) continue;

      double hsx1 = shower1->get_hsx();
      double hsy1 = shower1->get_hsy();
//      double hsdr = (hsx1-hsx0)*(hsx1-hsx0)+(hsy1-hsy0)*(hsy1-hsy0);
//      hsdr = sqrt(hsdr);
//      if(hsdr < 0.01) continue;
      ratio = shower1->get_esum()/shower1->get_mpcECorr();
      rms_x = shower1->get_rms_hsx();
      rms_y = shower1->get_rms_hsy();
      rms = sqrt(rms_x*rms_x+rms_y*rms_y);
      clus_dist = shower1->get_ClosestMPCClusterDistance();
      if(ratio < 0.015 || ratio > 0.9 || rms < 0.0015 || rms > 0.004|| clus_dist > 0.003) continue;

      double norm0 = sqrt(hsx0*hsx0+hsy0*hsy0+1*1);
      double norm1 = sqrt(hsx1*hsx1+hsy1*hsy1+1*1);
      double prdct = hsx0*hsx1+hsy0*hsy1+1*1;
      double cs = prdct/(norm0*norm1);
      double mass = sqrt(2*e0*e1*(1-cs));
//      cout <<"mass: "<<mass<<endl;
      hPi0_mass[arm0]->Fill(mass);
    }
   
  }
 
  return EVENT_OK;
}

int mMpcExPi0Reco::event_display(PHCompositeNode* topNode){
  fkinWrapper* fkin = findNode::getClass<fkinWrapper>(topNode, "fkin");
  primaryWrapper* primary = findNode::getClass<primaryWrapper>(topNode, "primary");
  
  
  if(!fkin || !primary){
    cout <<PHWHERE<<"No fkinWrapper !!!"<<endl;
    return EVENT_OK;
  }
  size_t fkinrows = fkin->RowCount();
  size_t primrows = fkin->RowCount();
  //it looks phi and theta is in degree
  double tot_e0 = 0;
  double theta0 = 0;
  double phi0 = 0;
  double tot_e1 = 0;
  double theta1 = 0;
  double phi1 = 0;
  bool isfirst = true;
  bool isPi0 = false;

  
  //no decay photon is recorded in primary
//  for(size_t iprim = 0;iprim < primrows;iprim++){
//    int idpart = primary->get_idpart(iprim);
//    cout <<idpart<<endl;
//    if(iprim > 3) break;
//  }
  
  for(size_t ifkin = 0;ifkin < fkinrows;ifkin++){
    int idpart = fkin->get_idpart(ifkin);
    int idparent = fkin->get_idparent(ifkin);
//    cout <<idpart<<" "<<idparent<<endl;
    //get the photon
    if(idpart != 1) continue;
    //pi0 
    if(idparent == 7){
//      cout <<"find pi0 !"<<endl;
      if(isfirst){
        tot_e0 = fkin->get_ptot(ifkin);
	theta0 = fkin->get_pthet(ifkin)*3.14159/180.;
	phi0 = fkin->get_pphi(ifkin)*3.14159/180.;
//	cout <<"first photon: "<<tot_e0<<" "<<theta0<<" "<<phi0<<endl;
	isfirst = false;
      }
      else{
        tot_e1 = fkin->get_ptot(ifkin);
	theta1 = fkin->get_pthet(ifkin)*3.14159/180.;
	phi1 = fkin->get_pphi(ifkin)*3.14159/180.;
	isPi0 = true;
//        cout <<"second photon: "<<tot_e1<<" "<<theta1<<" "<<phi1<<endl;
	break;
      }
    }
  }
 
  if(tot_e0+tot_e1 < 1 || tot_e0+tot_e1 > 25) return EVENT_OK;
  cout <<"photon 1 "<<tot_e0<<" photon 2 "<<tot_e1<<endl;
  for(int i = 0;i < 2; i++){
    if(grammyl[i]) delete grammyl[i];
    if(grammyh[i]) delete grammyh[i];
    char name[50];
    sprintf(name,"Arm %d for low q",i);
    grammyl[i] = new Exogram(name,name,900,-24,24,900,-24,24,8,-0.5,7.5);
    sprintf(name,"Arm %d for high q",i);
    grammyh[i] = new Exogram(name,name,900,-24,24,900,-24,24,8,-0.5,7.5);
  }

  int Nhits = _mpcex_hit_container->size();
  cout << "Number of Hits: "<<Nhits<<endl;
  for(int i = 0;i < Nhits;i++){
    TMpcExHit* hit = _mpcex_hit_container->getHit(i);
    double high_q = hit->high();
    double low_q = hit->low();
    int key = hit->key();
    int arm = hit->arm();
    if(hit->isGoodHighHit()){
      if(high_q > 5) grammyh[arm]->FillEx(key,high_q);
    }
    if(hit->isGoodLowHit()){
      if(low_q > 5) grammyl[arm]->FillEx(key,low_q);
    }
  }

  for(int i = 0;i < 2;i++){
    if(hmpc_gridxy[i]) delete hmpc_gridxy[i];
    char name[100];
    sprintf(name,"hmpc_gridxy_arm%d",i);
    hmpc_gridxy[i] =  new TH2D(name,name,600,-24,24,600,-24,24);
  }

 _mpc_map = MpcMap::instance();
  int NMpcTowers = _mpc_tower_container->size();
cout<<"NMpcTowers: "<<NMpcTowers<<endl;
  for(int itower = 0; itower < NMpcTowers;itower++){
    mpcTowerContent* ctwr = _mpc_tower_container->getTower(itower);
    int tow_ch = ctwr->get_ch();
    int arm = _mpc_map->getArm(tow_ch);
    double e_tower = ctwr->get_energy();
    if(e_tower < 0.) continue;
    cout<<"e_tower "<<e_tower<<endl;
    double x = _mpc_map->getX(tow_ch);
    double y = _mpc_map->getY(tow_ch);
    int binx0 = hmpc_gridxy[arm]->GetXaxis()->FindBin(x-0.9);
    int binx1 = hmpc_gridxy[arm]->GetXaxis()->FindBin(x+0.9);
    int biny0 = hmpc_gridxy[arm]->GetYaxis()->FindBin(y-0.9);
    int biny1 = hmpc_gridxy[arm]->GetYaxis()->FindBin(y+0.9);
    for(int i = binx0;i <=binx1;i++){
      for(int j = biny0;j <=biny1;j++){
        hmpc_gridxy[arm]->SetBinContent(i,j,e_tower);
      }
    }
  }

  TCanvas* c_mpc0 = new TCanvas("c_mpc0","c_mpc0",1500,800);
  c_mpc0->Divide(2,1);
  c_mpc0->cd(1);
  hmpc_gridxy[0]->Draw("colz");
  c_mpc0->cd(2);
  grammyh[0]->Project3D("yx")->DrawCopy("colz");
  TCanvas* c_mpc1 = new TCanvas("c_mpc1","c_mpc1",1500,800);
  c_mpc1->Divide(2,1);
  c_mpc1->cd(1);
  hmpc_gridxy[1]->Draw("colz");
  c_mpc1->cd(2);
  grammyh[1]->Project3D("yx")->DrawCopy("colz");
  
  return EVENT_OK;
}

struct MpcShower{
  int arm;
  double mpc_e;
  double mpcex_e;
  double hx;
  double hy;
};

int mMpcExPi0Reco::pi0_mass_simV2_study(PHCompositeNode* topNode){
  fkinWrapper* fkin = findNode::getClass<fkinWrapper>(topNode, "fkin");
  primaryWrapper* primary = findNode::getClass<primaryWrapper>(topNode, "primary");
  
  if(!fkin || !primary){
    cout <<PHWHERE<<"No fkinWrapper !!!"<<endl;
    return EVENT_OK;
  }
  size_t fkinrows = fkin->RowCount();
  size_t primrows = fkin->RowCount();
  //it looks phi and theta is in degree
  double tot_e0 = 0;
  double theta0 = 0;
  double phi0 = 0;
  double tot_e1 = 0;
  double theta1 = 0;
  double phi1 = 0;
  bool isfirst = true;
  bool isPi0 = false;

  
  for(size_t ifkin = 0;ifkin < fkinrows;ifkin++){
    int idpart = fkin->get_idpart(ifkin);
    int idparent = fkin->get_idparent(ifkin);
//    cout <<idpart<<" "<<idparent<<endl;
    //get the photon
    if(idpart != 1) continue;
    //pi0 
    if(idparent == 7){
//      cout <<"find pi0 !"<<endl;
      if(isfirst){
        tot_e0 = fkin->get_ptot(ifkin);
	theta0 = fkin->get_pthet(ifkin)*3.14159/180.;
	phi0 = fkin->get_pphi(ifkin)*3.14159/180.;
//	cout <<"first photon: "<<tot_e0<<" "<<theta0<<" "<<phi0<<endl;
	isfirst = false;
      }
      else{
        tot_e1 = fkin->get_ptot(ifkin);
	theta1 = fkin->get_pthet(ifkin)*3.14159/180.;
	phi1 = fkin->get_pphi(ifkin)*3.14159/180.;
	isPi0 = true;
//        cout <<"second photon: "<<tot_e1<<" "<<theta1<<" "<<phi1<<endl;
	break;
      }
    }
  }

  

  if(tot_e0+tot_e1>25 || tot_e0+tot_e1<0.5) return EVENT_OK;
  
  int Ntowers = _mpc_tower_container->size();
//cout <<"Ntowers: "<<Ntowers<<endl;
  int Nshowers = _mpcex_shower_container->size();
  set<int>used_showers;
  _mpc_map = MpcMap::instance();
  vector<MpcShower> mpc_shower_list;
  for(int i = 0;i < Ntowers;i++){
    mpcTowerContent* tower = _mpc_tower_container->getTower(i);
    int tower_ch = tower->get_ch();
//cout<<"tower ch "<<tower_ch<<endl;
    int gridx = _mpc_map->getGridX(tower_ch);
    int gridy = _mpc_map->getGridY(tower_ch);
    int arm = _mpc_map->getArm(tower_ch);
    double tower_x = _mpc_map->getX(tower_ch);
    double tower_y = _mpc_map->getY(tower_ch);
//cout<<"tower x "<<tower_x<<" tower y "<<tower_y<<endl;
    double e = tower->get_energy();
//cout<<"tower e: "<<e<<endl;
    if(e < 0) continue;
    double E3x3 = 0;
    double grx = 0;
    double gry = 0;
    double grx2 = 0;
    double gry2 = 0;
    bool ispeak = true;
    int tower_count = 0;
    for(int ix = gridx-1;ix <= gridx+1;ix++){
      for(int iy = gridy-1;iy <=gridy+1;iy++){
        int temp_ch = _mpc_map->getFeeCh(ix,iy,arm);
//cout<<"ix: "<<ix<<" iy: "<<iy<<" temp ch: "<<temp_ch<<endl;
	if(temp_ch > 576 || temp_ch < 0) continue;
	int temp_index = _mpc_tower_container->findTower(temp_ch);
//cout<<"temp index: "<<temp_index<<endl;
	if(temp_index < 0) continue;
	mpcTowerContent* temp_twr = _mpc_tower_container->getTower(temp_index);
        double temp_e = temp_twr->get_energy();
//cout<<"temp e: "<<temp_e<<endl;
        if(temp_e < 0) continue;
	if(temp_e > e){
          ispeak = false;
	  break;
	}
	tower_count++;
	E3x3 += temp_e;
	double temp_x = _mpc_map->getX(temp_ch);
	double temp_y = _mpc_map->getY(temp_ch);
//cout<<"temp x "<<temp_x<<endl;
//cout<<"temp y "<<temp_y<<endl;
	grx += temp_x*temp_e;
	gry += temp_y*temp_e;
	grx2 += temp_x*temp_x*temp_e;
	gry2 += temp_y*temp_y*temp_e;
      }
      if(!ispeak) break;
    }
    if(!ispeak) continue;
//    cout <<"find peak tower "<<" E3x3: "<<E3x3<<endl;
    if(e/E3x3 > 0.95 || tower_count < 3) continue;
    grx = grx/E3x3;
    gry = gry/E3x3;
//cout<<"grx: "<<grx<<" gry: "<<gry<<endl;
    grx2 = grx2/E3x3;
    gry2 = gry2/E3x3;
    double lz = 220.9 - _vertex;
    if(arm == 0) lz = -220.9 - _vertex;
    double ghx = grx/lz;
    double ghy = gry/lz;
//cout<<"ghx: "<<ghx<<" ghy: "<<ghy<<endl; 
    double grms_x = sqrt(grx2-grx*grx);
    if(grx2 <= grx*grx || grms_x < 0.00001) grms_x = 1.1/sqrt(12.);
    double grms_y = sqrt(gry2-gry*gry);
    if(gry2 <= gry*gry || grms_y < 0.00001) grms_y = 1.1/sqrt(12.);
    
    double gh_rms_x = fabs(grms_x/lz);
    double gh_rms_y = fabs(grms_y/lz);
//cout<<"gh_rms_x: "<<gh_rms_x<<" gh_rms_y: "<<gh_rms_y<<endl;
    double mpcex_hx = 0;
    double mpcex_hy = 0;
    double mpcex_esum = 0;
    int shower_count = 0;
    for(int j = 0;j < Nshowers;j++){
//      if(used_showers.find(j) != used_showers.end()) continue;
      TMpcExShower* shower = _mpcex_shower_container->getShower(j);
      int sh_arm = shower->get_arm();
      if(sh_arm != arm ) continue;
      double sh_hsx = shower->get_hsx();
      double sh_hsy = shower->get_hsy();
//cout <<"shower hsx: "<<sh_hsx<<"shower hsy: "<<sh_hsy<<endl;
      double dhx = sh_hsx - ghx;
      double dhy = sh_hsy - ghy;
//cout <<"dhx: "<<dhx<<" dhy: "<<dhy<<endl;
      if(fabs(dhx) < 3*gh_rms_x && fabs(dhy) < 3*gh_rms_y){
        double sh_e = shower->get_esum();
	if(sh_e < 0) continue;
	mpcex_hx += sh_hsx*sh_e;
	mpcex_hy += sh_hsy*sh_e;
	mpcex_esum += sh_e;
        used_showers.insert(j);
	shower_count++;
      }
    }
    if(shower_count == 0){
//      cout <<"no shower found for the tower !!!"<<endl;
      continue;
    }
    mpcex_hx = mpcex_hx/mpcex_esum;
    mpcex_hy = mpcex_hy/mpcex_esum;
    MpcShower mpc_shower;
    mpc_shower.mpc_e = E3x3;
    mpc_shower.mpcex_e = mpcex_esum;
    mpc_shower.hx = mpcex_hx;
    mpc_shower.hy = mpcex_hy;
    mpc_shower.arm = arm;
    mpc_shower_list.push_back(mpc_shower);
//cout <<"shower count: "<<shower_count<<endl;
//cout <<" mpc_e: "<<E3x3<<" mpcex_e: "<<mpcex_esum<<" total_e "<<E3x3+mpcex_esum
//     <<" hx: "<<mpcex_hx<<" hy: "<<mpcex_hy<<endl
//     <<" x: "<<mpcex_hx*lz
//     <<" y: "<<mpcex_hy*lz
//     <<endl;

  }//i tower
  if(mpc_shower_list.size() < 2) return EVENT_OK;
  int Nsize = mpc_shower_list.size();
  for(int i = 0;i < Nsize;i++){
    MpcShower mpc_shower0 = mpc_shower_list[i];
    double e0 = mpc_shower0.mpc_e + mpc_shower0.mpcex_e;
    double hsx0 = mpc_shower0.hx;
    double hsy0 = mpc_shower0.hy;
    int arm0 = mpc_shower0.arm;
    for(int j = i+1;j < Nsize;j++){
      MpcShower mpc_shower1 = mpc_shower_list[j];
      double e1 = mpc_shower1.mpc_e+mpc_shower1.mpcex_e;
      if(e0+e1> 25 || e0+e1<1) continue;
      double hsx1 = mpc_shower1.hx;
      double hsy1 = mpc_shower1.hy;
      int arm1 = mpc_shower1.arm;
      if(arm0 != arm1) continue;
      double norm0 = sqrt(hsx0*hsx0+hsy0*hsy0+1*1);
      double norm1 = sqrt(hsx1*hsx1+hsy1*hsy1+1*1);
      double prdct = hsx0*hsx1+hsy0*hsy1+1*1;
      double cs = prdct/(norm0*norm1);
      double mass = sqrt(2*e0*e1*(1-cs));
//      cout <<"mass: "<<mass<<endl;
      hPi0_mass[0]->Fill(mass);
    }
  }

  return EVENT_OK;
}
