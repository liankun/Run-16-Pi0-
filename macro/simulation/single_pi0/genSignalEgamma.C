TCanvas *Canvas;
TH1F *hEta;
TH1F *hE;

void genSignalEgamma(Int_t fNEvents=1000,
		     Float_t fLowE=5, Float_t fHighE = 5,
		     Float_t fLowEta=3.5, Float_t fHighEta=3.5)
{
  hEta = new TH1F("hEta","",100,-5,5);
  hE = new TH1F("hPt","",100,0,50);
  gRandom->SetSeed(gSystem->GetPid());
  Float_t fMass = 0;
  ofstream fileOut("oscar.input");
  fileOut << "# OSC1999A" << endl;
  fileOut << "# final_id_p_x" << endl;
  fileOut << "# exodus event generator in single particle mode" << endl;
  fileOut << "#" << endl;
  for(Int_t i=0 ; i<fNEvents ; i++) {
    Float_t fEnergy = gRandom->Rndm()*(fHighE-fLowE)+fLowE;
    Float_t fPhi = (gRandom->Rndm()-0.5)*2*TMath::Pi();
    Float_t fEta = pow(-1,i)*(gRandom->Rndm()*(fHighEta-fLowEta)+fLowEta);

    Float_t fP  = TMath::Power(fEnergy*fEnergy-fMass*fMass,0.5);
    Float_t fTheta = atan(exp(fabs(fEta)))*2.;
    Float_t fPz = TMath::Sign(fP*TMath::Cos(fTheta),fEta);
    Float_t fPt = TMath::Power(fP*fP-fPz*fPz,0.5);
    
    //Create the needed variables: px, py, pz and Energy
    Float_t fPx = fPt*TMath::Cos(fPhi);
    Float_t fPy = fPt*TMath::Sin(fPhi);
    //cout << fEnergy << " " << TMath::Power(fPz*fPz+fPt*fPt+fMass*fMass,0.5) << endl;
    hEta->Fill(fEta);
    hPt->Fill(fPt);
    fileOut << "0 1" << endl;
    fileOut << "0 22 0 " << fPx << " " << fPy << " " << fPz << " "
	    << fEnergy << " 0 0 0 0 0" << endl;
    fileOut << "0 0" << endl;}
  return;
  gROOT->SetStyle("Plain");
  Canvas = new TCanvas("Canvas","",0,0,800,400);
  Canvas->Divide(2,1);
  Canvas->cd(1); hEta->Draw();
  Canvas->cd(2); hPt->Draw();
} 
