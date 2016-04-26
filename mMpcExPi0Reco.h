#ifndef __MPCEXPI0RECO_HH__
#define __MPCEXPI0RECO_HH__

#ifndef __CINT__
#include <SubsysReco.h>
#endif

class SubsysReco;
class PHCompositeNode; 
class MpcMap;
class mpcClusterContainer;
class mpcTowerContainer;
class TMpcExHitContainer;
class TMpcExShowerContainer;
class TMpcExLShowerContainer;
class TMpcExHit;
class TH1D;
class TH2D;
class TH3D;
class Exogram;

class mMpcExPi0Reco:public SubsysReco{
  public:
    mMpcExPi0Reco(const char* name = "MPCEXPI0RECO");
    virtual int Init(PHCompositeNode*);
    virtual int InitRun(PHCompositeNode*);
    virtual int process_event(PHCompositeNode*);
    virtual ~mMpcExPi0Reco();
    virtual int End(PHCompositeNode*);

    int pi0_mass_sim_init();
    int pi0_mass_simV2_study(PHCompositeNode*);
    int pi0_mass_sim(PHCompositeNode*); //John's shower module ,single Pi0
    int pi0_mass_sim_PP(PHCompositeNode*);//PP pythia
    int pi0_mass_data(PHCompositeNode*);
    int event_display(PHCompositeNode*);
    private:
      void set_interface_ptrs(PHCompositeNode*);
      TMpcExHitContainer* _mpcex_hit_container;
      TMpcExShowerContainer* _mpcex_shower_container;
      TMpcExLShowerContainer* _mpcex_lshower_container;
      mpcClusterContainer* _mpc_cluster_container;
      mpcTowerContainer* _mpc_tower_container;
      MpcMap* _mpc_map;


      TH1D* hPi0_mass[2];
      TH2D* hmass_dr[2];
      TH2D* hmass_angle[2];//mass vs open angle
      TH3D* hAsy_tote_angle;
      TH2D* he_reco_real0[2];//first photon
      TH2D* he_reco_real1[2];//second photon
      TH2D* hnomatch_e_ratio[2];//energy ratio between MpcEx and MPC
      TH2D* hnomatch_rms[2];
      TH2D* hnomatch_cluster[2];//the closest distance to cluster
      TH2D* hnomatch_fired_layers;
      TH2D* hnomatch_Nhits;
      TH2D* htote_reco_real;
      TH2D* hangle_reco_real;
      TH2D* htote_clus;//merge the shower if they have same cluster
      TH2D* hangle_clus;//merge the shower if they have the same cluster
      TH1D* hmass_merge;
 
      
      //event display
      Exogram* grammyh[2];
      Exogram* grammyl[2];
      TH2D* hmpc_gridxy[2];

    protected:
      double _vertex;
};

#endif 
