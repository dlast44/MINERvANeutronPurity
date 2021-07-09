//File: EventLoop.cxx
//Info: This is a script to run a loop over all events in a single nTuple file and perform some plotting. Will eventually exist as the basis for the loops over events in analysis.
//
//Usage: EventLoop.cxx <MasterAnaDev_NTuple_list/single_file> <0=MC/1=PC> <0=tracker/1=targets/2=both> <0=trueSignalOnly/1=trueBackgroundOnly/2=all> <output_directory> <tag_for_naming_files> optional: <n_event g.t. 0 if you want constraint otherwise it'll do all> <1="Dan's",anything else default> <PC non-muon EnergyCut>
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

//C++ includes
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <bitset>
#include <time.h>
#include <sys/stat.h>

//ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TMath.h"

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"

#include "syst/CVUniverse.h"
#include "obj/NeutCands.h"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

using namespace std;

double targetBoundary = 5850.0;
bitset<4> goodBlob{"1111"};

bool PassesFVCuts(CVUniverse& univ, int region){
  vector<double> vtx = univ.GetVtx();
  double side = 850.0*2.0/sqrt(3.0);
  if (region == 0){
    if (vtx[2] < targetBoundary || vtx[2] > 8422.0 ) return false;
  }
  else if (region == 1){
    if (vtx[2] > targetBoundary || vtx[2] < 4500.0 ) return false;
  }
  else if (vtx[2] < 4500.0 || vtx[2] > 8422.0) return false;
  if (fabs(vtx[0]) > 850.00) return false;
  double slope = 1.0/sqrt(3.0);
  
  if(fabs(vtx[1]) < side - slope*fabs(vtx[0]))return true;
  return false;
}

bool PassesCleanCCAntiNuCuts(CVUniverse& univ, int isPC, double ECut=10000.0){
  int MINOSMatch=0;
  if (isPC){
    //    MINOSMatch=univ.GetIsMinosMatchTrackOLD();
    MINOSMatch=1;
    TLorentzVector prim_part(univ.GetFSPartPx().at(0),univ.GetFSPartPy().at(0),univ.GetFSPartPz().at(0),univ.GetFSPartE().at(0));
    if ((prim_part.E()-prim_part.M()) > ECut) return false;
  }
  else{
    MINOSMatch=univ.GetHasInteractionVertex();
      //MINOSMatch=univ.GetIsMinosMatchTrack();
  }
  return
    (univ.GetNDeadDiscriminatorsUpstreamMuon() < 2) &&
    (univ.GetNuHelicity() == 2) &&
    (MINOSMatch == 1) &&
    (TMath::RadToDeg()*univ.GetThetamu() < 20.0) &&
    (univ.GetPmu() < 20000.0 && univ.GetPmu() > 1500.0);
}

//Should Code the more general anti-nu CCQE cuts at some point... but this is a focus for Tejin stuff...
bool PassesTejinRecoilCut(CVUniverse& univ, int isPC, int whichRecoil){
  if (isPC) return true;
  double Q2GeV = univ.GetQ2QEPickledGeV();
  double recoilEGeV = univ.GetRecoilEnergyGeV();
  if (whichRecoil == 1) recoilEGeV = univ.GetDANRecoilEnergyGeV();
  if (Q2GeV < 0.0 || recoilEGeV < 0.0) return false;
  else if (Q2GeV < 0.3) return (recoilEGeV < (0.04+0.43*Q2GeV));
  else if (Q2GeV < 1.4) return (recoilEGeV < (0.08+0.3*Q2GeV));
  else return (recoilEGeV < 0.5);
}

bool PassesTejinCCQECuts(CVUniverse& univ){
  //bool PassesRecoilECut = false;
  //recoil energy cut and Q2 calculation are unclear to me. Need investigate...
  vector<double> EMBlobInfo = univ.GetEMNBlobsTotalEnergyTotalNHits();
  return 
    (univ.GetNTracks() == 1) &&
    (EMBlobInfo.at(0) < 2) &&
    (EMBlobInfo.at(1) >= 10.0*EMBlobInfo.at(2)) &&
    !(univ.GetNImprovedMichel() > 0);//Should this be !=0 ??????????
}

int PassesTejinBlobCuts(CVUniverse& univ){
  NeutronCandidates::NeutCand leadingBlob=univ.GetCurrentLeadingNeutCand();
  int leading3D=leadingBlob.GetIs3D();
  double leadingEGeV = 0.001*leadingBlob.GetTotalE();
  TVector3 leadingPos=leadingBlob.GetBegPos();
  TVector3 leadingFP=leadingBlob.GetFlightPath();
  TVector3 muonMom(univ.GetMuon4V().Z(),univ.GetMuon4V().Y(),univ.GetMuon4V().Z());
  double Q2GeV = univ.GetQ2QEPickledGeV();
  double MnGeV = 0.939566; //Should check what others might use, but this will be close enough for now...
  if (leadingFP.Mag()==0 || muonMom.Mag()==0) return 0;
  else{
    if((leading3D==1) && (leadingFP.Angle(muonMom) > 0.261799388)){
      if (leadingPos.Z() >= targetBoundary) return 2;
      return 1;
    }
    return 0;
  }
}

bool PassesCuts(CVUniverse& univ, int isPC, int region, double ECut=10000.0){
  if (ECut <= 0.0){
    ECut = 10000.0;
  }
  return
    PassesFVCuts(univ, region) &&
    PassesCleanCCAntiNuCuts(univ, isPC, ECut) &&
    PassesTejinCCQECuts(univ);
}

bool IsTrueSignal(CVUniverse& univ){
  int current = univ.GetMCCurrent();
  int incoming = univ.GetMCIncoming();

  if (current != 1 || incoming != -14) return false;

  int genie_n_muons = 0;
  int genie_n_mesons = 0;
  int genie_n_heavy_baryons_plus_pi0s = 0;
  int genie_n_photons =0;
  int genie_n_protons = 0;

  vector<int> PDGs = univ.GetFSPartPDG();
  vector<double> Es = univ.GetFSPartE();

  for (unsigned int i=0; i<PDGs.size(); ++i){
    int pdg = PDGs.at(i);
    double energy = Es.at(i);
    double proton_E = 1058.272;
    if (abs(pdg) == 13) genie_n_muons++;
    else if ( pdg == 22  && energy > 10) genie_n_photons++;
    else if ( abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 323 || pdg == 111 || pdg == 130 || pdg == 310 || pdg == 311 || pdg == 313 ){
      genie_n_mesons++;
    }
    else if ( pdg == 3112 || pdg == 3122 || pdg == 3212 || pdg == 3222 || pdg == 4112 || pdg == 4122 || pdg == 4222 || pdg == 411 || pdg == 421 || pdg == 111){
      genie_n_heavy_baryons_plus_pi0s++;
    }
    else if ( pdg == 2212 && energy > proton_E ) genie_n_protons++;
  }

  return genie_n_muons == 1 && 
    genie_n_mesons == 0 && 
    genie_n_heavy_baryons_plus_pi0s == 0 && 
    genie_n_photons == 0 && 
    genie_n_protons == 0;
}

bool PathExists(string path){
  struct stat buffer;
  return (stat (path.c_str(), &buffer) == 0);
}

void Write1DHistsToFile(vector<PlotUtils::HistWrapper<CVUniverse>*> hists, CVUniverse* univ, TFile* file){
  for (auto hist : hists){
    hist->SyncCVHistos();
    hist->univHist(univ)->SetDirectory(file);
    hist->univHist(univ)->Write();
  }
}

int main(int argc, char* argv[]) {

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc < 7 || argc > 10) {
    cout << "Check usage..." << endl;
    return 2;
  }

  string playlist=string(argv[1]);
  int isPC=atoi(argv[2]);
  int region=atoi(argv[3]);
  int sample=atoi(argv[4]);
  string outDir=string(argv[5]);
  string tag=string(argv[6]);
  int nEntries=0;
  int whichRecoil=0;
  double PCECut=-1.0;

  if (argc >= 8){
    nEntries=atoi(argv[7]);
  }
  if (argc>=9){
    whichRecoil=atoi(argv[8]);
  }
  if (argc == 10){
    PCECut=atof(argv[9]);
  }

  if (PathExists(outDir)){
    cout << "Thank you for choosing a path for output files that exists." << endl;
  }
  else{
    cout << "Output directory doesn't exist. Exiting" << endl;
    return 3;
  }

  if (region < 0 || region > 2){
    cout << "Check usage for meaning of different regions." << endl;
    return 4;
  } 

  if (sample < 0 || sample > 2){
    cout << "Check usage for meaning of different samples." << endl;
    return 5;
  } 

  string txtExt = ".txt";
  string rootExt = ".root";
  string slash = "/";
  string token;
  string playlistStub = playlist;
  size_t pos=0;

  //cout << playlistStub << endl;
  while ((pos = playlistStub.find(slash)) != string::npos){
    //cout << playlistStub << endl;
    token = playlistStub.substr(0, pos);
    //cout << token << endl;
    playlistStub.erase(0, pos+slash.length());
  }
  //cout << playlistStub << endl;
  if ((pos=playlistStub.find(txtExt)) != string::npos){
    token = playlistStub.substr(0,pos);
  }
  else if ((pos=playlistStub.find(rootExt)) != string::npos){
    token = playlistStub.substr(0,pos);
  }
  else{
    cout << "input must either be .root, or .txt" << endl;
    return 3;
  }

  playlistStub=token;
  cout << "Input file name parsed to: " << playlistStub << endl;

  map<int,TString>regionNames={{0,"tracker"},{1,"nuke"},{2,"fullID"},};
  map<int,TString>sampleNames={{0,"Signal"},{1,"Background"},{2, "AllSelected"}};

  unordered_map<int,int> PDGbins;
  PDGbins[2112] = 2;
  PDGbins[2212] = 3;
  PDGbins[111] = 4;
  PDGbins[211] = 5;
  PDGbins[-211] = 6;
  PDGbins[22] = 7;
  PDGbins[11] = 8;
  PDGbins[-11] = 8;
  PDGbins[-13] = 9;
  PDGbins[13] = 9;

  PlotUtils::ChainWrapper* chain = makeChainWrapperPtr(playlist,"MasterAnaDev");
  
  CVUniverse* CV = new CVUniverse(chain);
  map< string, vector<CVUniverse*>> error_bands;
  error_bands[string("CV")].push_back(CV);

  map<int, TString> typeNames={{1,"QE"},{2,"RES"},{3,"DIS"},{8,"2p2h"},{0,"Other"}};

  //Blob Level Plots
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_primary_parent_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_primary_parent_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_primary_parent_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_primary_parent_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_length_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_length_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_length_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_length_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_avg_dEdx_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_avg_dEdx_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_avg_dEdx_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_avg_dEdx_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_blobE_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_blobE_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_blobE_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_blobE_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_dist_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_dist_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_dist_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_dist_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_Zdist_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_Zdist_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_Zdist_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_tracker_Zdist_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_primary_parent_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_primary_parent_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_primary_parent_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_primary_parent_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_length_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_length_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_length_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_length_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_avg_dEdx_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_avg_dEdx_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_avg_dEdx_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_avg_dEdx_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_blobE_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_blobE_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_blobE_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_blobE_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_dist_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_dist_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_dist_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_dist_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_Zdist_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_Zdist_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_Zdist_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_target_Zdist_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_primary_parent_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_primary_parent_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_primary_parent_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_primary_parent_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_length_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_length_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_length_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_length_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_avg_dEdx_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_avg_dEdx_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_avg_dEdx_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_avg_dEdx_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_blobE_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_blobE_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_blobE_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_blobE_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_dist_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_dist_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_dist_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_dist_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_Zdist_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_Zdist_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_Zdist_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_ALL_Zdist_Tejin_TrackerONLY;

  //Event Level Plots
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_leadBlob_passes_classifier_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_leadBlob_passes_classifier_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_leadBlob_passes_classifier_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_leadBlob_passes_classifier_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_n3DBlobs_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_n3DBlobs_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_n3DBlobs_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_n3DBlobs_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_nGoodBlobs_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_nGoodBlobs_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_nGoodBlobs_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_nGoodBlobs_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_nBlobs_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_nBlobs_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_nBlobs_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_nBlobs_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_AvgBlobEnergy_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_AvgBlobEnergy_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_AvgBlobEnergy_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_AvgBlobEnergy_Tejin_TrackerONLY;

  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_RecoilEnergyGeV_CCQE;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_RecoilEnergyGeV_Recoil;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_RecoilEnergyGeV_Tejin;
  map<int, PlotUtils::HistWrapper<CVUniverse>> map_hw_RecoilEnergyGeV_Tejin_TrackerONLY;

  for(auto type: typeNames){
    map_hw_tracker_primary_parent_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_primary_parent_CCQE_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE);;Blobs",10,0,10,error_bands);
    map_hw_tracker_primary_parent_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_primary_parent_Recoil_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil);;Blobs",10,0,10,error_bands);
    map_hw_tracker_primary_parent_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_primary_parent_Tejin_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil, Blob);;Blobs",10,0,10,error_bands);
    map_hw_tracker_primary_parent_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_primary_parent_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil, Blob, Tracker Blob);;Blobs",10,0,10,error_bands);

    map_hw_tracker_length_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_length_CCQE_"+type.second,"True "+type.second+" Blob Length (CCQE);Len. [mm];Blobs",50,0,500,error_bands);
    map_hw_tracker_length_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_length_Recoil_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil);Len. [mm];Blobs",50,0,500,error_bands);
    map_hw_tracker_length_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_length_Tejin_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil, Blob);Len. [mm];Blobs",50,0,500,error_bands);
    map_hw_tracker_length_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_length_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil, Blob, Tracker Blob);Len. [mm];Blobs",50,0,500,error_bands);

    map_hw_tracker_avg_dEdx_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_avg_dEdx_CCQE_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE);dE/dx [MeV/mm];Blobs",25,0,50,error_bands);
    map_hw_tracker_avg_dEdx_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_avg_dEdx_Recoil_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil);dE/dx [MeV/mm];Blobs",25,0,50,error_bands);
    map_hw_tracker_avg_dEdx_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_avg_dEdx_Tejin_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil, Blob);dE/dx [MeV/mm];Blobs",25,0,50,error_bands);
    map_hw_tracker_avg_dEdx_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_avg_dEdx_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil, Blob, Tracker Blob);dE/dx [MeV/mm];Blobs",25,0,50,error_bands);

    map_hw_tracker_blobE_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_blobE_CCQE_"+type.second,"True "+type.second+" Blob Energy (CCQE);E [MeV];Blobs",50,0,150,error_bands);
    map_hw_tracker_blobE_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_blobE_Recoil_"+type.second,"True "+type.second+" Blob Energy (CCQE, Recoil);E [MeV];Blobs",50,0,150,error_bands);
    map_hw_tracker_blobE_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_blobE_Tejin_"+type.second,"True "+type.second+" Blob Energy (CCQE, Recoil, Blob);E [MeV];Blobs",50,0,150,error_bands);
    map_hw_tracker_blobE_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_blobE_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Energy (CCQE, Recoil, Blob, Tracker Blob);E [MeV];Blobs",50,0,150,error_bands);

    map_hw_tracker_dist_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_dist_CCQE_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_tracker_dist_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_dist_Recoil_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_tracker_dist_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_dist_Tejin_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil, Blob);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_tracker_dist_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_dist_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil, Blob, Tracker Blob);Dist. [mm];Blobs",300,0,3000,error_bands);

    map_hw_tracker_Zdist_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_Zdist_CCQE_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_tracker_Zdist_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_Zdist_Recoil_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_tracker_Zdist_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_Zdist_Tejin_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil, Blob);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_tracker_Zdist_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_Zdist_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil, Blob, Tracker Blob);Dist. [mm];Blobs",300,0,3000,error_bands);

    map_hw_target_primary_parent_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_primary_parent_CCQE_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE);;Blobs",10,0,10,error_bands);
    map_hw_target_primary_parent_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_primary_parent_Recoil_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil);;Blobs",10,0,10,error_bands);
    map_hw_target_primary_parent_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_primary_parent_Tejin_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil, Blob);;Blobs",10,0,10,error_bands);
    map_hw_target_primary_parent_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_primary_parent_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil, Blob, Tracker Blob);;Blobs",10,0,10,error_bands);

    map_hw_target_length_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_length_CCQE_"+type.second,"True "+type.second+" Blob Length (CCQE);Len. [mm];Blobs",50,0,500,error_bands);
    map_hw_target_length_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_length_Recoil_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil);Len. [mm];Blobs",50,0,500,error_bands);
    map_hw_target_length_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_length_Tejin_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil, Blob);Len. [mm];Blobs",50,0,500,error_bands);
    map_hw_target_length_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_length_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil, Blob, Tracker Blob);Len. [mm];Blobs",50,0,500,error_bands);

    map_hw_target_avg_dEdx_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_avg_dEdx_CCQE_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE);dE/dx [MeV/mm];Blobs",25,0,50,error_bands);
    map_hw_target_avg_dEdx_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_avg_dEdx_Recoil_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil);dE/dx [MeV/mm];Blobs",25,0,50,error_bands);
    map_hw_target_avg_dEdx_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_avg_dEdx_Tejin_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil, Blob);dE/dx [MeV/mm];Blobs",25,0,50,error_bands);
    map_hw_target_avg_dEdx_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_avg_dEdx_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil, Blob, Tracker Blob);dE/dx [MeV/mm];Blobs",25,0,50,error_bands);

    map_hw_target_blobE_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_blobE_CCQE_"+type.second,"True "+type.second+" Blob Energy (CCQE);E [MeV];Blobs",50,0,150,error_bands);
    map_hw_target_blobE_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_blobE_Recoil_"+type.second,"True "+type.second+" Blob Energy (CCQE, Recoil);E [MeV];Blobs",50,0,150,error_bands);
    map_hw_target_blobE_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_blobE_Tejin_"+type.second,"True "+type.second+" Blob Energy (CCQE, Recoil, Blob);E [MeV];Blobs",50,0,150,error_bands);
    map_hw_target_blobE_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_blobE_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Energy (CCQE, Recoil, Blob, Tracker Blob);E [MeV];Blobs",50,0,150,error_bands);

    map_hw_target_dist_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_dist_CCQE_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_target_dist_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_dist_Recoil_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_target_dist_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_dist_Tejin_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil, Blob);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_target_dist_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_dist_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil, Blob, Tracker Blob);Dist. [mm];Blobs",300,0,3000,error_bands);

    map_hw_target_Zdist_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_Zdist_CCQE_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_target_Zdist_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_Zdist_Recoil_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_target_Zdist_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_Zdist_Tejin_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil, Blob);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_target_Zdist_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_Zdist_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil, Blob, Tracker Blob);Dist. [mm];Blobs",300,0,3000,error_bands);

    map_hw_ALL_primary_parent_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_primary_parent_CCQE_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE);;Blobs",10,0,10,error_bands);
    map_hw_ALL_primary_parent_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_primary_parent_Recoil_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil);;Blobs",10,0,10,error_bands);
    map_hw_ALL_primary_parent_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_primary_parent_Tejin_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil, Blob);;Blobs",10,0,10,error_bands);
    map_hw_ALL_primary_parent_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_primary_parent_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil, Blob, Tracker Blob);;Blobs",10,0,10,error_bands);

    map_hw_ALL_length_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_length_CCQE_"+type.second,"True "+type.second+" Blob Length (CCQE);Len. [mm];Blobs",50,0,500,error_bands);
    map_hw_ALL_length_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_length_Recoil_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil);Len. [mm];Blobs",50,0,500,error_bands);
    map_hw_ALL_length_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_length_Tejin_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil, Blob);Len. [mm];Blobs",50,0,500,error_bands);
    map_hw_ALL_length_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_length_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil, Blob, Tracker Blob);Len. [mm];Blobs",50,0,500,error_bands);

    map_hw_ALL_avg_dEdx_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_avg_dEdx_CCQE_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE);dE/dx [MeV/mm];Blobs",25,0,50,error_bands);
    map_hw_ALL_avg_dEdx_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_avg_dEdx_Recoil_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil);dE/dx [MeV/mm];Blobs",25,0,50,error_bands);
    map_hw_ALL_avg_dEdx_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_avg_dEdx_Tejin_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil, Blob);dE/dx [MeV/mm];Blobs",25,0,50,error_bands);
    map_hw_ALL_avg_dEdx_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_avg_dEdx_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil, Blob, Tracker Blob);dE/dx [MeV/mm];Blobs",25,0,50,error_bands);

    map_hw_ALL_blobE_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_blobE_CCQE_"+type.second,"True "+type.second+" Blob Energy (CCQE);E [MeV];Blobs",50,0,150,error_bands);
    map_hw_ALL_blobE_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_blobE_Recoil_"+type.second,"True "+type.second+" Blob Energy (CCQE, Recoil);E [MeV];Blobs",50,0,150,error_bands);
    map_hw_ALL_blobE_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_blobE_Tejin_"+type.second,"True "+type.second+" Blob Energy (CCQE, Recoil, Blob);E [MeV];Blobs",50,0,150,error_bands);
    map_hw_ALL_blobE_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_blobE_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Energy (CCQE, Recoil, Blob, Tracker Blob);E [MeV];Blobs",50,0,150,error_bands);

    map_hw_ALL_dist_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_dist_CCQE_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_ALL_dist_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_dist_Recoil_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_ALL_dist_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_dist_Tejin_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil, Blob);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_ALL_dist_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_dist_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil, Blob, Tracker Blob);Dist. [mm];Blobs",300,0,3000,error_bands);

    map_hw_ALL_Zdist_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_Zdist_CCQE_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_ALL_Zdist_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_Zdist_Recoil_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_ALL_Zdist_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_Zdist_Tejin_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil, Blob);Dist. [mm];Blobs",300,0,3000,error_bands);
    map_hw_ALL_Zdist_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_Zdist_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil, Blob, Tracker Blob);Dist. [mm];Blobs",300,0,3000,error_bands);

    map_hw_leadBlob_passes_classifier_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_leadBlob_passes_classifier_CCQE_"+type.second,"True "+type.second+" Leading Blob Passes (CCQE);;Events",2,0,2,error_bands);
    map_hw_leadBlob_passes_classifier_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_leadBlob_passes_classifier_Recoil_"+type.second,"True "+type.second+" Leading Blob Passes (CCQE, Recoil);;Events",2,0,2,error_bands);
    map_hw_leadBlob_passes_classifier_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_leadBlob_passes_classifier_Tejin_"+type.second,"True "+type.second+" Leading Blob Passes (CCQE, Recoil, Blob);;Events",2,0,2,error_bands);
    map_hw_leadBlob_passes_classifier_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_leadBlob_passes_classifier_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Leading Blob Passes (CCQE, Recoil, Blob, Tracker Blob);;Events",2,0,2,error_bands);

    map_hw_n3DBlobs_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_n3DBlobs_CCQE_"+type.second,"True "+type.second+" No. 3D Blobs (CCQE);No.;Events",10,0,10,error_bands);
    map_hw_n3DBlobs_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_n3DBlobs_Recoil_"+type.second,"True "+type.second+" No. 3D Blobs (CCQE, Recoil);No.;Events",10,0,10,error_bands);
    map_hw_n3DBlobs_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_n3DBlobs_Tejin_"+type.second,"True "+type.second+" No. 3D Blobs (CCQE, Recoil, Blob);No.;Events",10,0,10,error_bands);
    map_hw_n3DBlobs_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_n3DBlobs_Tejin_TrackerONLY_"+type.second,"True "+type.second+" No. 3D Blobs (CCQE, Recoil, Blob, Tracker Blob);No.;Events",10,0,10,error_bands);

    map_hw_nGoodBlobs_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nGoodBlobs_CCQE_"+type.second,"True "+type.second+" No. Blobs Which Pass (CCQE);No.;Events",10,0,10,error_bands);
    map_hw_nGoodBlobs_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nGoodBlobs_Recoil_"+type.second,"True "+type.second+" No. Blobs Which Pass (CCQE, Recoil);No.;Events",10,0,10,error_bands);
    map_hw_nGoodBlobs_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nGoodBlobs_Tejin_"+type.second,"True "+type.second+" No. Blobs Which Pass (CCQE, Recoil, Blob);No.;Events",10,0,10,error_bands);
    map_hw_nGoodBlobs_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nGoodBlobs_Tejin_TrackerONLY_"+type.second,"True "+type.second+" No. Blobs Which Pass (CCQE, Recoil, Blob, Tracker Blob);No.;Events",10,0,10,error_bands);

    map_hw_nBlobs_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nBlobs_CCQE_"+type.second,"True "+type.second+" No. of Blobs (CCQE);No.;Events",100,0,100,error_bands);
    map_hw_nBlobs_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nBlobs_Recoil_"+type.second,"True "+type.second+" No. of Blobs (CCQE, Recoil);No.;Events",100,0,100,error_bands);
    map_hw_nBlobs_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nBlobs_Tejin_"+type.second,"True "+type.second+" No. of Blobs (CCQE, Recoil, Blob);No.;Events",100,0,100,error_bands);
    map_hw_nBlobs_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nBlobs_Tejin_TrackerONLY_"+type.second,"True "+type.second+" No. of Blobs (CCQE, Recoil, Blob, Tracker Blob);No.;Events",100,0,100,error_bands);

    map_hw_AvgBlobEnergy_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_AvgBlobEnergy_CCQE_"+type.second,"True "+type.second+" Avg. Blob Energy (CCQE);Avg. E [MeV];Events",50,0,50,error_bands);
    map_hw_AvgBlobEnergy_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_AvgBlobEnergy_Recoil_"+type.second,"True "+type.second+" Avg. Blob Energy (CCQE, Recoil);Avg. E [MeV];Events",50,0,50,error_bands);
    map_hw_AvgBlobEnergy_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_AvgBlobEnergy_Tejin_"+type.second,"True "+type.second+" Avg. Blob Energy (CCQE, Recoil, Blob);Avg. E [MeV];Events",50,0,50,error_bands);
    map_hw_AvgBlobEnergy_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_AvgBlobEnergy_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Avg. Blob Energy (CCQE, Recoil, Blob, Tracker Blob);Avg. E [MeV];Events",50,0,50,error_bands);

    map_hw_RecoilEnergyGeV_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_RecoilEnergyGeV_CCQE_"+type.second,"True "+type.second+" Recoil Energy (CCQE);RecoilE [GeV];Events",50,0,1.5,error_bands);
    map_hw_RecoilEnergyGeV_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_RecoilEnergyGeV_Recoil_"+type.second,"True "+type.second+" Recoil Energy (CCQE, Recoil);RecoilE [GeV];Events",50,0,1.5,error_bands);
    map_hw_RecoilEnergyGeV_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_RecoilEnergyGeV_Tejin_"+type.second,"True "+type.second+" Recoil Energy (CCQE, Recoil, Blob);RecoilE [GeV];Events",50,0,1.5,error_bands);
    map_hw_RecoilEnergyGeV_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_RecoilEnergyGeV_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Recoil Energy (CCQE, Recoil, Blob, Tracker Blob);RecoilE [GeV];Events",50,0,1.5,error_bands);
  }

  vector<PlotUtils::HistWrapper<CVUniverse>*> histsALL ={
    &map_hw_tracker_primary_parent_CCQE[1],&map_hw_tracker_primary_parent_CCQE[2],&map_hw_tracker_primary_parent_CCQE[3],&map_hw_tracker_primary_parent_CCQE[8],&map_hw_tracker_primary_parent_CCQE[0],
    &map_hw_tracker_primary_parent_Recoil[1],&map_hw_tracker_primary_parent_Recoil[2],&map_hw_tracker_primary_parent_Recoil[3],&map_hw_tracker_primary_parent_Recoil[8],&map_hw_tracker_primary_parent_Recoil[0],
    &map_hw_tracker_primary_parent_Tejin[1],&map_hw_tracker_primary_parent_Tejin[2],&map_hw_tracker_primary_parent_Tejin[3],&map_hw_tracker_primary_parent_Tejin[8],&map_hw_tracker_primary_parent_Tejin[0],
    &map_hw_tracker_primary_parent_Tejin_TrackerONLY[1],&map_hw_tracker_primary_parent_Tejin_TrackerONLY[2],&map_hw_tracker_primary_parent_Tejin_TrackerONLY[3],&map_hw_tracker_primary_parent_Tejin_TrackerONLY[8],&map_hw_tracker_primary_parent_Tejin_TrackerONLY[0],

    &map_hw_tracker_length_CCQE[1],&map_hw_tracker_length_CCQE[2],&map_hw_tracker_length_CCQE[3],&map_hw_tracker_length_CCQE[8],&map_hw_tracker_length_CCQE[0],
    &map_hw_tracker_length_Recoil[1],&map_hw_tracker_length_Recoil[2],&map_hw_tracker_length_Recoil[3],&map_hw_tracker_length_Recoil[8],&map_hw_tracker_length_Recoil[0],
    &map_hw_tracker_length_Tejin[1],&map_hw_tracker_length_Tejin[2],&map_hw_tracker_length_Tejin[3],&map_hw_tracker_length_Tejin[8],&map_hw_tracker_length_Tejin[0],
    &map_hw_tracker_length_Tejin_TrackerONLY[1],&map_hw_tracker_length_Tejin_TrackerONLY[2],&map_hw_tracker_length_Tejin_TrackerONLY[3],&map_hw_tracker_length_Tejin_TrackerONLY[8],&map_hw_tracker_length_Tejin_TrackerONLY[0],

    &map_hw_tracker_avg_dEdx_CCQE[1],&map_hw_tracker_avg_dEdx_CCQE[2],&map_hw_tracker_avg_dEdx_CCQE[3],&map_hw_tracker_avg_dEdx_CCQE[8],&map_hw_tracker_avg_dEdx_CCQE[0],
    &map_hw_tracker_avg_dEdx_Recoil[1],&map_hw_tracker_avg_dEdx_Recoil[2],&map_hw_tracker_avg_dEdx_Recoil[3],&map_hw_tracker_avg_dEdx_Recoil[8],&map_hw_tracker_avg_dEdx_Recoil[0],
    &map_hw_tracker_avg_dEdx_Tejin[1],&map_hw_tracker_avg_dEdx_Tejin[2],&map_hw_tracker_avg_dEdx_Tejin[3],&map_hw_tracker_avg_dEdx_Tejin[8],&map_hw_tracker_avg_dEdx_Tejin[0],
    &map_hw_tracker_avg_dEdx_Tejin_TrackerONLY[1],&map_hw_tracker_avg_dEdx_Tejin_TrackerONLY[2],&map_hw_tracker_avg_dEdx_Tejin_TrackerONLY[3],&map_hw_tracker_avg_dEdx_Tejin_TrackerONLY[8],&map_hw_tracker_avg_dEdx_Tejin_TrackerONLY[0],

    &map_hw_tracker_blobE_CCQE[1],&map_hw_tracker_blobE_CCQE[2],&map_hw_tracker_blobE_CCQE[3],&map_hw_tracker_blobE_CCQE[8],&map_hw_tracker_blobE_CCQE[0],
    &map_hw_tracker_blobE_Recoil[1],&map_hw_tracker_blobE_Recoil[2],&map_hw_tracker_blobE_Recoil[3],&map_hw_tracker_blobE_Recoil[8],&map_hw_tracker_blobE_Recoil[0],
    &map_hw_tracker_blobE_Tejin[1],&map_hw_tracker_blobE_Tejin[2],&map_hw_tracker_blobE_Tejin[3],&map_hw_tracker_blobE_Tejin[8],&map_hw_tracker_blobE_Tejin[0],
    &map_hw_tracker_blobE_Tejin_TrackerONLY[1],&map_hw_tracker_blobE_Tejin_TrackerONLY[2],&map_hw_tracker_blobE_Tejin_TrackerONLY[3],&map_hw_tracker_blobE_Tejin_TrackerONLY[8],&map_hw_tracker_blobE_Tejin_TrackerONLY[0],

    &map_hw_tracker_dist_CCQE[1],&map_hw_tracker_dist_CCQE[2],&map_hw_tracker_dist_CCQE[3],&map_hw_tracker_dist_CCQE[8],&map_hw_tracker_dist_CCQE[0],
    &map_hw_tracker_dist_Recoil[1],&map_hw_tracker_dist_Recoil[2],&map_hw_tracker_dist_Recoil[3],&map_hw_tracker_dist_Recoil[8],&map_hw_tracker_dist_Recoil[0],
    &map_hw_tracker_dist_Tejin[1],&map_hw_tracker_dist_Tejin[2],&map_hw_tracker_dist_Tejin[3],&map_hw_tracker_dist_Tejin[8],&map_hw_tracker_dist_Tejin[0],
    &map_hw_tracker_dist_Tejin_TrackerONLY[1],&map_hw_tracker_dist_Tejin_TrackerONLY[2],&map_hw_tracker_dist_Tejin_TrackerONLY[3],&map_hw_tracker_dist_Tejin_TrackerONLY[8],&map_hw_tracker_dist_Tejin_TrackerONLY[0],

    &map_hw_tracker_Zdist_CCQE[1],&map_hw_tracker_Zdist_CCQE[2],&map_hw_tracker_Zdist_CCQE[3],&map_hw_tracker_Zdist_CCQE[8],&map_hw_tracker_Zdist_CCQE[0],
    &map_hw_tracker_Zdist_Recoil[1],&map_hw_tracker_Zdist_Recoil[2],&map_hw_tracker_Zdist_Recoil[3],&map_hw_tracker_Zdist_Recoil[8],&map_hw_tracker_Zdist_Recoil[0],
    &map_hw_tracker_Zdist_Tejin[1],&map_hw_tracker_Zdist_Tejin[2],&map_hw_tracker_Zdist_Tejin[3],&map_hw_tracker_Zdist_Tejin[8],&map_hw_tracker_Zdist_Tejin[0],
    &map_hw_tracker_Zdist_Tejin_TrackerONLY[1],&map_hw_tracker_Zdist_Tejin_TrackerONLY[2],&map_hw_tracker_Zdist_Tejin_TrackerONLY[3],&map_hw_tracker_Zdist_Tejin_TrackerONLY[8],&map_hw_tracker_Zdist_Tejin_TrackerONLY[0],

    &map_hw_target_primary_parent_CCQE[1],&map_hw_target_primary_parent_CCQE[2],&map_hw_target_primary_parent_CCQE[3],&map_hw_target_primary_parent_CCQE[8],&map_hw_target_primary_parent_CCQE[0],
    &map_hw_target_primary_parent_Recoil[1],&map_hw_target_primary_parent_Recoil[2],&map_hw_target_primary_parent_Recoil[3],&map_hw_target_primary_parent_Recoil[8],&map_hw_target_primary_parent_Recoil[0],
    &map_hw_target_primary_parent_Tejin[1],&map_hw_target_primary_parent_Tejin[2],&map_hw_target_primary_parent_Tejin[3],&map_hw_target_primary_parent_Tejin[8],&map_hw_target_primary_parent_Tejin[0],
    &map_hw_target_primary_parent_Tejin_TrackerONLY[1],&map_hw_target_primary_parent_Tejin_TrackerONLY[2],&map_hw_target_primary_parent_Tejin_TrackerONLY[3],&map_hw_target_primary_parent_Tejin_TrackerONLY[8],&map_hw_target_primary_parent_Tejin_TrackerONLY[0],

    &map_hw_target_length_CCQE[1],&map_hw_target_length_CCQE[2],&map_hw_target_length_CCQE[3],&map_hw_target_length_CCQE[8],&map_hw_target_length_CCQE[0],
    &map_hw_target_length_Recoil[1],&map_hw_target_length_Recoil[2],&map_hw_target_length_Recoil[3],&map_hw_target_length_Recoil[8],&map_hw_target_length_Recoil[0],
    &map_hw_target_length_Tejin[1],&map_hw_target_length_Tejin[2],&map_hw_target_length_Tejin[3],&map_hw_target_length_Tejin[8],&map_hw_target_length_Tejin[0],
    &map_hw_target_length_Tejin_TrackerONLY[1],&map_hw_target_length_Tejin_TrackerONLY[2],&map_hw_target_length_Tejin_TrackerONLY[3],&map_hw_target_length_Tejin_TrackerONLY[8],&map_hw_target_length_Tejin_TrackerONLY[0],

    &map_hw_target_avg_dEdx_CCQE[1],&map_hw_target_avg_dEdx_CCQE[2],&map_hw_target_avg_dEdx_CCQE[3],&map_hw_target_avg_dEdx_CCQE[8],&map_hw_target_avg_dEdx_CCQE[0],
    &map_hw_target_avg_dEdx_Recoil[1],&map_hw_target_avg_dEdx_Recoil[2],&map_hw_target_avg_dEdx_Recoil[3],&map_hw_target_avg_dEdx_Recoil[8],&map_hw_target_avg_dEdx_Recoil[0],
    &map_hw_target_avg_dEdx_Tejin[1],&map_hw_target_avg_dEdx_Tejin[2],&map_hw_target_avg_dEdx_Tejin[3],&map_hw_target_avg_dEdx_Tejin[8],&map_hw_target_avg_dEdx_Tejin[0],
    &map_hw_target_avg_dEdx_Tejin_TrackerONLY[1],&map_hw_target_avg_dEdx_Tejin_TrackerONLY[2],&map_hw_target_avg_dEdx_Tejin_TrackerONLY[3],&map_hw_target_avg_dEdx_Tejin_TrackerONLY[8],&map_hw_target_avg_dEdx_Tejin_TrackerONLY[0],

    &map_hw_target_blobE_CCQE[1],&map_hw_target_blobE_CCQE[2],&map_hw_target_blobE_CCQE[3],&map_hw_target_blobE_CCQE[8],&map_hw_target_blobE_CCQE[0],
    &map_hw_target_blobE_Recoil[1],&map_hw_target_blobE_Recoil[2],&map_hw_target_blobE_Recoil[3],&map_hw_target_blobE_Recoil[8],&map_hw_target_blobE_Recoil[0],
    &map_hw_target_blobE_Tejin[1],&map_hw_target_blobE_Tejin[2],&map_hw_target_blobE_Tejin[3],&map_hw_target_blobE_Tejin[8],&map_hw_target_blobE_Tejin[0],
    &map_hw_target_blobE_Tejin_TrackerONLY[1],&map_hw_target_blobE_Tejin_TrackerONLY[2],&map_hw_target_blobE_Tejin_TrackerONLY[3],&map_hw_target_blobE_Tejin_TrackerONLY[8],&map_hw_target_blobE_Tejin_TrackerONLY[0],

    &map_hw_target_dist_CCQE[1],&map_hw_target_dist_CCQE[2],&map_hw_target_dist_CCQE[3],&map_hw_target_dist_CCQE[8],&map_hw_target_dist_CCQE[0],
    &map_hw_target_dist_Recoil[1],&map_hw_target_dist_Recoil[2],&map_hw_target_dist_Recoil[3],&map_hw_target_dist_Recoil[8],&map_hw_target_dist_Recoil[0],
    &map_hw_target_dist_Tejin[1],&map_hw_target_dist_Tejin[2],&map_hw_target_dist_Tejin[3],&map_hw_target_dist_Tejin[8],&map_hw_target_dist_Tejin[0],
    &map_hw_target_dist_Tejin_TrackerONLY[1],&map_hw_target_dist_Tejin_TrackerONLY[2],&map_hw_target_dist_Tejin_TrackerONLY[3],&map_hw_target_dist_Tejin_TrackerONLY[8],&map_hw_target_dist_Tejin_TrackerONLY[0],

    &map_hw_target_Zdist_CCQE[1],&map_hw_target_Zdist_CCQE[2],&map_hw_target_Zdist_CCQE[3],&map_hw_target_Zdist_CCQE[8],&map_hw_target_Zdist_CCQE[0],
    &map_hw_target_Zdist_Recoil[1],&map_hw_target_Zdist_Recoil[2],&map_hw_target_Zdist_Recoil[3],&map_hw_target_Zdist_Recoil[8],&map_hw_target_Zdist_Recoil[0],
    &map_hw_target_Zdist_Tejin[1],&map_hw_target_Zdist_Tejin[2],&map_hw_target_Zdist_Tejin[3],&map_hw_target_Zdist_Tejin[8],&map_hw_target_Zdist_Tejin[0],
    &map_hw_target_Zdist_Tejin_TrackerONLY[1],&map_hw_target_Zdist_Tejin_TrackerONLY[2],&map_hw_target_Zdist_Tejin_TrackerONLY[3],&map_hw_target_Zdist_Tejin_TrackerONLY[8],&map_hw_target_Zdist_Tejin_TrackerONLY[0],

    &map_hw_ALL_primary_parent_CCQE[1],&map_hw_ALL_primary_parent_CCQE[2],&map_hw_ALL_primary_parent_CCQE[3],&map_hw_ALL_primary_parent_CCQE[8],&map_hw_ALL_primary_parent_CCQE[0],
    &map_hw_ALL_primary_parent_Recoil[1],&map_hw_ALL_primary_parent_Recoil[2],&map_hw_ALL_primary_parent_Recoil[3],&map_hw_ALL_primary_parent_Recoil[8],&map_hw_ALL_primary_parent_Recoil[0],
    &map_hw_ALL_primary_parent_Tejin[1],&map_hw_ALL_primary_parent_Tejin[2],&map_hw_ALL_primary_parent_Tejin[3],&map_hw_ALL_primary_parent_Tejin[8],&map_hw_ALL_primary_parent_Tejin[0],
    &map_hw_ALL_primary_parent_Tejin_TrackerONLY[1],&map_hw_ALL_primary_parent_Tejin_TrackerONLY[2],&map_hw_ALL_primary_parent_Tejin_TrackerONLY[3],&map_hw_ALL_primary_parent_Tejin_TrackerONLY[8],&map_hw_ALL_primary_parent_Tejin_TrackerONLY[0],

    &map_hw_ALL_length_CCQE[1],&map_hw_ALL_length_CCQE[2],&map_hw_ALL_length_CCQE[3],&map_hw_ALL_length_CCQE[8],&map_hw_ALL_length_CCQE[0],
    &map_hw_ALL_length_Recoil[1],&map_hw_ALL_length_Recoil[2],&map_hw_ALL_length_Recoil[3],&map_hw_ALL_length_Recoil[8],&map_hw_ALL_length_Recoil[0],
    &map_hw_ALL_length_Tejin[1],&map_hw_ALL_length_Tejin[2],&map_hw_ALL_length_Tejin[3],&map_hw_ALL_length_Tejin[8],&map_hw_ALL_length_Tejin[0],
    &map_hw_ALL_length_Tejin_TrackerONLY[1],&map_hw_ALL_length_Tejin_TrackerONLY[2],&map_hw_ALL_length_Tejin_TrackerONLY[3],&map_hw_ALL_length_Tejin_TrackerONLY[8],&map_hw_ALL_length_Tejin_TrackerONLY[0],

    &map_hw_ALL_avg_dEdx_CCQE[1],&map_hw_ALL_avg_dEdx_CCQE[2],&map_hw_ALL_avg_dEdx_CCQE[3],&map_hw_ALL_avg_dEdx_CCQE[8],&map_hw_ALL_avg_dEdx_CCQE[0],
    &map_hw_ALL_avg_dEdx_Recoil[1],&map_hw_ALL_avg_dEdx_Recoil[2],&map_hw_ALL_avg_dEdx_Recoil[3],&map_hw_ALL_avg_dEdx_Recoil[8],&map_hw_ALL_avg_dEdx_Recoil[0],
    &map_hw_ALL_avg_dEdx_Tejin[1],&map_hw_ALL_avg_dEdx_Tejin[2],&map_hw_ALL_avg_dEdx_Tejin[3],&map_hw_ALL_avg_dEdx_Tejin[8],&map_hw_ALL_avg_dEdx_Tejin[0],
    &map_hw_ALL_avg_dEdx_Tejin_TrackerONLY[1],&map_hw_ALL_avg_dEdx_Tejin_TrackerONLY[2],&map_hw_ALL_avg_dEdx_Tejin_TrackerONLY[3],&map_hw_ALL_avg_dEdx_Tejin_TrackerONLY[8],&map_hw_ALL_avg_dEdx_Tejin_TrackerONLY[0],

    &map_hw_ALL_blobE_CCQE[1],&map_hw_ALL_blobE_CCQE[2],&map_hw_ALL_blobE_CCQE[3],&map_hw_ALL_blobE_CCQE[8],&map_hw_ALL_blobE_CCQE[0],
    &map_hw_ALL_blobE_Recoil[1],&map_hw_ALL_blobE_Recoil[2],&map_hw_ALL_blobE_Recoil[3],&map_hw_ALL_blobE_Recoil[8],&map_hw_ALL_blobE_Recoil[0],
    &map_hw_ALL_blobE_Tejin[1],&map_hw_ALL_blobE_Tejin[2],&map_hw_ALL_blobE_Tejin[3],&map_hw_ALL_blobE_Tejin[8],&map_hw_ALL_blobE_Tejin[0],
    &map_hw_ALL_blobE_Tejin_TrackerONLY[1],&map_hw_ALL_blobE_Tejin_TrackerONLY[2],&map_hw_ALL_blobE_Tejin_TrackerONLY[3],&map_hw_ALL_blobE_Tejin_TrackerONLY[8],&map_hw_ALL_blobE_Tejin_TrackerONLY[0],

    &map_hw_ALL_dist_CCQE[1],&map_hw_ALL_dist_CCQE[2],&map_hw_ALL_dist_CCQE[3],&map_hw_ALL_dist_CCQE[8],&map_hw_ALL_dist_CCQE[0],
    &map_hw_ALL_dist_Recoil[1],&map_hw_ALL_dist_Recoil[2],&map_hw_ALL_dist_Recoil[3],&map_hw_ALL_dist_Recoil[8],&map_hw_ALL_dist_Recoil[0],
    &map_hw_ALL_dist_Tejin[1],&map_hw_ALL_dist_Tejin[2],&map_hw_ALL_dist_Tejin[3],&map_hw_ALL_dist_Tejin[8],&map_hw_ALL_dist_Tejin[0],
    &map_hw_ALL_dist_Tejin_TrackerONLY[1],&map_hw_ALL_dist_Tejin_TrackerONLY[2],&map_hw_ALL_dist_Tejin_TrackerONLY[3],&map_hw_ALL_dist_Tejin_TrackerONLY[8],&map_hw_ALL_dist_Tejin_TrackerONLY[0],

    &map_hw_ALL_Zdist_CCQE[1],&map_hw_ALL_Zdist_CCQE[2],&map_hw_ALL_Zdist_CCQE[3],&map_hw_ALL_Zdist_CCQE[8],&map_hw_ALL_Zdist_CCQE[0],
    &map_hw_ALL_Zdist_Recoil[1],&map_hw_ALL_Zdist_Recoil[2],&map_hw_ALL_Zdist_Recoil[3],&map_hw_ALL_Zdist_Recoil[8],&map_hw_ALL_Zdist_Recoil[0],
    &map_hw_ALL_Zdist_Tejin[1],&map_hw_ALL_Zdist_Tejin[2],&map_hw_ALL_Zdist_Tejin[3],&map_hw_ALL_Zdist_Tejin[8],&map_hw_ALL_Zdist_Tejin[0],
    &map_hw_ALL_Zdist_Tejin_TrackerONLY[1],&map_hw_ALL_Zdist_Tejin_TrackerONLY[2],&map_hw_ALL_Zdist_Tejin_TrackerONLY[3],&map_hw_ALL_Zdist_Tejin_TrackerONLY[8],&map_hw_ALL_Zdist_Tejin_TrackerONLY[0],

    &map_hw_leadBlob_passes_classifier_CCQE[1],&map_hw_leadBlob_passes_classifier_CCQE[2],&map_hw_leadBlob_passes_classifier_CCQE[3],&map_hw_leadBlob_passes_classifier_CCQE[8],&map_hw_leadBlob_passes_classifier_CCQE[0],
    &map_hw_leadBlob_passes_classifier_Recoil[1],&map_hw_leadBlob_passes_classifier_Recoil[2],&map_hw_leadBlob_passes_classifier_Recoil[3],&map_hw_leadBlob_passes_classifier_Recoil[8],&map_hw_leadBlob_passes_classifier_Recoil[0],
    &map_hw_leadBlob_passes_classifier_Tejin[1],&map_hw_leadBlob_passes_classifier_Tejin[2],&map_hw_leadBlob_passes_classifier_Tejin[3],&map_hw_leadBlob_passes_classifier_Tejin[8],&map_hw_leadBlob_passes_classifier_Tejin[0],
    &map_hw_leadBlob_passes_classifier_Tejin_TrackerONLY[1],&map_hw_leadBlob_passes_classifier_Tejin_TrackerONLY[2],&map_hw_leadBlob_passes_classifier_Tejin_TrackerONLY[3],&map_hw_leadBlob_passes_classifier_Tejin_TrackerONLY[8],&map_hw_leadBlob_passes_classifier_Tejin_TrackerONLY[0],

    &map_hw_n3DBlobs_CCQE[1],&map_hw_n3DBlobs_CCQE[2],&map_hw_n3DBlobs_CCQE[3],&map_hw_n3DBlobs_CCQE[8],&map_hw_n3DBlobs_CCQE[0],
    &map_hw_n3DBlobs_Recoil[1],&map_hw_n3DBlobs_Recoil[2],&map_hw_n3DBlobs_Recoil[3],&map_hw_n3DBlobs_Recoil[8],&map_hw_n3DBlobs_Recoil[0],
    &map_hw_n3DBlobs_Tejin[1],&map_hw_n3DBlobs_Tejin[2],&map_hw_n3DBlobs_Tejin[3],&map_hw_n3DBlobs_Tejin[8],&map_hw_n3DBlobs_Tejin[0],
    &map_hw_n3DBlobs_Tejin_TrackerONLY[1],&map_hw_n3DBlobs_Tejin_TrackerONLY[2],&map_hw_n3DBlobs_Tejin_TrackerONLY[3],&map_hw_n3DBlobs_Tejin_TrackerONLY[8],&map_hw_n3DBlobs_Tejin_TrackerONLY[0],

    &map_hw_nGoodBlobs_CCQE[1],&map_hw_nGoodBlobs_CCQE[2],&map_hw_nGoodBlobs_CCQE[3],&map_hw_nGoodBlobs_CCQE[8],&map_hw_nGoodBlobs_CCQE[0],
    &map_hw_nGoodBlobs_Recoil[1],&map_hw_nGoodBlobs_Recoil[2],&map_hw_nGoodBlobs_Recoil[3],&map_hw_nGoodBlobs_Recoil[8],&map_hw_nGoodBlobs_Recoil[0],
    &map_hw_nGoodBlobs_Tejin[1],&map_hw_nGoodBlobs_Tejin[2],&map_hw_nGoodBlobs_Tejin[3],&map_hw_nGoodBlobs_Tejin[8],&map_hw_nGoodBlobs_Tejin[0],
    &map_hw_nGoodBlobs_Tejin_TrackerONLY[1],&map_hw_nGoodBlobs_Tejin_TrackerONLY[2],&map_hw_nGoodBlobs_Tejin_TrackerONLY[3],&map_hw_nGoodBlobs_Tejin_TrackerONLY[8],&map_hw_nGoodBlobs_Tejin_TrackerONLY[0],

    &map_hw_nBlobs_CCQE[1],&map_hw_nBlobs_CCQE[2],&map_hw_nBlobs_CCQE[3],&map_hw_nBlobs_CCQE[8],&map_hw_nBlobs_CCQE[0],
    &map_hw_nBlobs_Recoil[1],&map_hw_nBlobs_Recoil[2],&map_hw_nBlobs_Recoil[3],&map_hw_nBlobs_Recoil[8],&map_hw_nBlobs_Recoil[0],
    &map_hw_nBlobs_Tejin[1],&map_hw_nBlobs_Tejin[2],&map_hw_nBlobs_Tejin[3],&map_hw_nBlobs_Tejin[8],&map_hw_nBlobs_Tejin[0],
    &map_hw_nBlobs_Tejin_TrackerONLY[1],&map_hw_nBlobs_Tejin_TrackerONLY[2],&map_hw_nBlobs_Tejin_TrackerONLY[3],&map_hw_nBlobs_Tejin_TrackerONLY[8],&map_hw_nBlobs_Tejin_TrackerONLY[0],

    &map_hw_AvgBlobEnergy_CCQE[1],&map_hw_AvgBlobEnergy_CCQE[2],&map_hw_AvgBlobEnergy_CCQE[3],&map_hw_AvgBlobEnergy_CCQE[8],&map_hw_AvgBlobEnergy_CCQE[0],
    &map_hw_AvgBlobEnergy_Recoil[1],&map_hw_AvgBlobEnergy_Recoil[2],&map_hw_AvgBlobEnergy_Recoil[3],&map_hw_AvgBlobEnergy_Recoil[8],&map_hw_AvgBlobEnergy_Recoil[0],
    &map_hw_AvgBlobEnergy_Tejin[1],&map_hw_AvgBlobEnergy_Tejin[2],&map_hw_AvgBlobEnergy_Tejin[3],&map_hw_AvgBlobEnergy_Tejin[8],&map_hw_AvgBlobEnergy_Tejin[0],
    &map_hw_AvgBlobEnergy_Tejin_TrackerONLY[1],&map_hw_AvgBlobEnergy_Tejin_TrackerONLY[2],&map_hw_AvgBlobEnergy_Tejin_TrackerONLY[3],&map_hw_AvgBlobEnergy_Tejin_TrackerONLY[8],&map_hw_AvgBlobEnergy_Tejin_TrackerONLY[0],

    &map_hw_RecoilEnergyGeV_CCQE[1],&map_hw_RecoilEnergyGeV_CCQE[2],&map_hw_RecoilEnergyGeV_CCQE[3],&map_hw_RecoilEnergyGeV_CCQE[8],&map_hw_RecoilEnergyGeV_CCQE[0],
    &map_hw_RecoilEnergyGeV_Recoil[1],&map_hw_RecoilEnergyGeV_Recoil[2],&map_hw_RecoilEnergyGeV_Recoil[3],&map_hw_RecoilEnergyGeV_Recoil[8],&map_hw_RecoilEnergyGeV_Recoil[0],
    &map_hw_RecoilEnergyGeV_Tejin[1],&map_hw_RecoilEnergyGeV_Tejin[2],&map_hw_RecoilEnergyGeV_Tejin[3],&map_hw_RecoilEnergyGeV_Tejin[8],&map_hw_RecoilEnergyGeV_Tejin[0],
    &map_hw_RecoilEnergyGeV_Tejin_TrackerONLY[1],&map_hw_RecoilEnergyGeV_Tejin_TrackerONLY[2],&map_hw_RecoilEnergyGeV_Tejin_TrackerONLY[3],&map_hw_RecoilEnergyGeV_Tejin_TrackerONLY[8],&map_hw_RecoilEnergyGeV_Tejin_TrackerONLY[0]
  };

  if(nEntries <= 0) nEntries = chain->GetEntries();
  cout << "Processing " << nEntries << " events." << endl;
  int n3DBlobs=0;
  int nGoodBlobs=0;
  double blobESum=0.0;
  for (int i=0; i<nEntries;++i){
    if (i%(nEntries/100)==0) cout << (100*i)/nEntries << "% finished." << endl;
    //if (i%(10000)==0) cout << i << " entries finished." << endl;
    for (auto band : error_bands){
      vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes){
	n3DBlobs=0;
	nGoodBlobs=0;
	blobESum=0.0;
	universe->SetEntry(i);
	universe->UpdateNeutCands();
	int nBlobs = universe->GetNNeutCands();
	double recoilEnergy = universe->GetRecoilEnergyGeV();
	if (whichRecoil==1) recoilEnergy = universe->GetDANRecoilEnergyGeV();
	if (sample == 0 && !IsTrueSignal(*universe)) continue;
	else if (sample == 1 && IsTrueSignal(*universe)) continue;
	//Passes CCQE Cuts that matche Tejin's selection
	if (PassesCuts(*universe, isPC, region, PCECut)){
	  
	  //int nFSPart = universe->GetNFSPart();
	  int intType = universe->GetInteractionType();
	  NeutronCandidates::NeutCand leadBlob = universe->GetCurrentLeadingNeutCand();
	  bool leadBlobPasses = (leadBlob.GetClassifier()==goodBlob);
	  if (intType > 8){
	    intType=0;
	  }
	  else if (intType < 1){
	    intType=0;
	  }
	  else if (intType > 3 && intType < 8){
	    intType=0;
	  }

	  //Passes Tejin Recoil and Blob
	  if (PassesTejinRecoilCut(*universe, isPC, whichRecoil)){
	    
	    int TejinBlobValue = PassesTejinBlobCuts(*universe);
	    //Passes Tejin Recoil and Blob
	    if (TejinBlobValue){
	      NeutronCandidates::NeutCands cands = universe->GetCurrentNeutCands();	      

	      for (auto& cand: cands.GetCandidates()){
		if (cand.second.GetIs3D()==1)++n3DBlobs;
		if (cand.second.GetClassifier()==goodBlob)++nGoodBlobs;

		int PID = cand.second.GetMCPID();
		int TopPID = cand.second.GetTopMCPID();
		int PTrackID = cand.second.GetMCParentTrackID();

		TVector3 FP = cand.second.GetFlightPath();
		double length = cand.second.GetDirection().Mag();
		double dEdx = -1.0;
		double blobE = cand.second.GetTotalE();
		if (length > 0.0) dEdx = blobE/length;
		double vtxDist = FP.Mag();
		double vtxZDist = abs(FP.Z());

		blobESum += blobE;
		/*
		if (PTrackID > nFSPart){
		  //Add something like this to learn how often this happened? ++nMultiIntBlobs;
		  continue;
		  }*/
		double candZ = cand.second.GetBegPos().Z();
		if (PTrackID==0 && !isPC){
		  //Additional Requirement of the Chosen Blob being in the tracker only.

		  if (TejinBlobValue==2){
		    map_hw_ALL_primary_parent_Tejin_TrackerONLY[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_ALL_length_Tejin_TrackerONLY[intType].univHist(universe)->Fill(length);
		    map_hw_ALL_avg_dEdx_Tejin_TrackerONLY[intType].univHist(universe)->Fill(dEdx);
		    map_hw_ALL_blobE_Tejin_TrackerONLY[intType].univHist(universe)->Fill(blobE);
		    map_hw_ALL_dist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_ALL_Zdist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxZDist);

		    if (candZ > targetBoundary){
		      map_hw_tracker_primary_parent_Tejin_TrackerONLY[intType].univHist(universe)->Fill(PDGbins[PID]);
		      map_hw_tracker_length_Tejin_TrackerONLY[intType].univHist(universe)->Fill(length);
		      map_hw_tracker_avg_dEdx_Tejin_TrackerONLY[intType].univHist(universe)->Fill(dEdx);
		      map_hw_tracker_blobE_Tejin_TrackerONLY[intType].univHist(universe)->Fill(blobE);
		      map_hw_tracker_dist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxDist);
		      map_hw_tracker_Zdist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxZDist);
		    }

		    else {
		      map_hw_target_primary_parent_Tejin_TrackerONLY[intType].univHist(universe)->Fill(PDGbins[PID]);
		      map_hw_target_length_Tejin_TrackerONLY[intType].univHist(universe)->Fill(length);
		      map_hw_target_avg_dEdx_Tejin_TrackerONLY[intType].univHist(universe)->Fill(dEdx);
		      map_hw_target_blobE_Tejin_TrackerONLY[intType].univHist(universe)->Fill(blobE);
		      map_hw_target_dist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxDist);
		      map_hw_target_Zdist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxZDist);
		    }
		  }

		  map_hw_ALL_primary_parent_Tejin[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_ALL_length_Tejin[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_Tejin[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_blobE_Tejin[intType].univHist(universe)->Fill(blobE);
		  map_hw_ALL_dist_Tejin[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_Tejin[intType].univHist(universe)->Fill(vtxZDist);
		  map_hw_ALL_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_ALL_length_Recoil[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_blobE_Recoil[intType].univHist(universe)->Fill(blobE);
		  map_hw_ALL_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		  map_hw_ALL_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_ALL_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		  map_hw_ALL_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);

		  if (candZ > targetBoundary){
		    map_hw_tracker_primary_parent_Tejin[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_tracker_length_Tejin[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_Tejin[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_blobE_Tejin[intType].univHist(universe)->Fill(blobE);
		    map_hw_tracker_dist_Tejin[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_Tejin[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_tracker_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_tracker_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_blobE_Recoil[intType].univHist(universe)->Fill(blobE);
		    map_hw_tracker_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_tracker_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_tracker_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		    map_hw_tracker_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }

		  else{
		    map_hw_target_primary_parent_Tejin[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_target_length_Tejin[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_Tejin[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_blobE_Tejin[intType].univHist(universe)->Fill(blobE);
		    map_hw_target_dist_Tejin[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_Tejin[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_target_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_target_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_blobE_Recoil[intType].univHist(universe)->Fill(blobE);
		    map_hw_target_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_target_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_target_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		    map_hw_target_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }
		}

		else{
		  //Additional Requirement of the Chosen Blob being in the tracker only.
		  if (TejinBlobValue==2){
		    map_hw_ALL_primary_parent_Tejin_TrackerONLY[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_ALL_length_Tejin_TrackerONLY[intType].univHist(universe)->Fill(length);
		    map_hw_ALL_avg_dEdx_Tejin_TrackerONLY[intType].univHist(universe)->Fill(dEdx);
		    map_hw_ALL_blobE_Tejin_TrackerONLY[intType].univHist(universe)->Fill(blobE);
		    map_hw_ALL_dist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_ALL_Zdist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxZDist);

		    if (candZ > targetBoundary){
		      map_hw_tracker_primary_parent_Tejin_TrackerONLY[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		      map_hw_tracker_length_Tejin_TrackerONLY[intType].univHist(universe)->Fill(length);
		      map_hw_tracker_avg_dEdx_Tejin_TrackerONLY[intType].univHist(universe)->Fill(dEdx);
		      map_hw_tracker_blobE_Tejin_TrackerONLY[intType].univHist(universe)->Fill(blobE);
		      map_hw_tracker_dist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxDist);
		      map_hw_tracker_Zdist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxZDist);
		    }		  

		    else {
		      map_hw_target_primary_parent_Tejin_TrackerONLY[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		      map_hw_target_length_Tejin_TrackerONLY[intType].univHist(universe)->Fill(length);
		      map_hw_target_avg_dEdx_Tejin_TrackerONLY[intType].univHist(universe)->Fill(dEdx);
		      map_hw_target_blobE_Tejin_TrackerONLY[intType].univHist(universe)->Fill(blobE);
		      map_hw_target_dist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxDist);
		      map_hw_target_Zdist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxZDist);
		    }
		  }

		  map_hw_ALL_primary_parent_Tejin[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_ALL_length_Tejin[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_Tejin[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_blobE_Tejin[intType].univHist(universe)->Fill(blobE);
		  map_hw_ALL_dist_Tejin[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_Tejin[intType].univHist(universe)->Fill(vtxZDist);
		  map_hw_ALL_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_ALL_length_Recoil[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_blobE_Recoil[intType].univHist(universe)->Fill(blobE);
		  map_hw_ALL_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		  map_hw_ALL_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_ALL_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		  map_hw_ALL_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);

		  if (candZ > targetBoundary){
		    map_hw_tracker_primary_parent_Tejin[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_tracker_length_Tejin[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_Tejin[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_blobE_Tejin[intType].univHist(universe)->Fill(blobE);
		    map_hw_tracker_dist_Tejin[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_Tejin[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_tracker_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_tracker_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_blobE_Recoil[intType].univHist(universe)->Fill(blobE);
		    map_hw_tracker_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_tracker_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_tracker_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		    map_hw_tracker_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }

		  else{
		    map_hw_target_primary_parent_Tejin[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_target_length_Tejin[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_Tejin[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_blobE_Tejin[intType].univHist(universe)->Fill(blobE);
		    map_hw_target_dist_Tejin[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_Tejin[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_target_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_target_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_blobE_Recoil[intType].univHist(universe)->Fill(blobE);
		    map_hw_target_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_target_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_target_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		    map_hw_target_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }
		}
	      }

	      //Event Level Plots for Passing Tejin Blob Cuts
	      if (leadBlobPasses) map_hw_leadBlob_passes_classifier_Tejin[intType].univHist(universe)->Fill(1);
	      else map_hw_leadBlob_passes_classifier_Tejin[intType].univHist(universe)->Fill(0);
	      map_hw_n3DBlobs_Tejin[intType].univHist(universe)->Fill(n3DBlobs);
	      map_hw_nGoodBlobs_Tejin[intType].univHist(universe)->Fill(nGoodBlobs);
	      map_hw_nBlobs_Tejin[intType].univHist(universe)->Fill(nBlobs);
	      map_hw_AvgBlobEnergy_Tejin[intType].univHist(universe)->Fill(blobESum/((double)(nBlobs)));
	      map_hw_RecoilEnergyGeV_Tejin[intType].univHist(universe)->Fill(recoilEnergy);

	      if (TejinBlobValue==2){
		if (leadBlobPasses) map_hw_leadBlob_passes_classifier_Tejin_TrackerONLY[intType].univHist(universe)->Fill(1);
		else map_hw_leadBlob_passes_classifier_Tejin_TrackerONLY[intType].univHist(universe)->Fill(0);
		map_hw_n3DBlobs_Tejin_TrackerONLY[intType].univHist(universe)->Fill(n3DBlobs);
		map_hw_nGoodBlobs_Tejin_TrackerONLY[intType].univHist(universe)->Fill(nGoodBlobs);
		map_hw_nBlobs_Tejin_TrackerONLY[intType].univHist(universe)->Fill(nBlobs);
		map_hw_AvgBlobEnergy_Tejin_TrackerONLY[intType].univHist(universe)->Fill(blobESum/((double)(nBlobs)));
		map_hw_RecoilEnergyGeV_Tejin_TrackerONLY[intType].univHist(universe)->Fill(recoilEnergy);
	      }
	    }
	  
	    //Passes Tejin Recoil Not Blob
	    else {
	      NeutronCandidates::NeutCands cands = universe->GetCurrentNeutCands();
	      for (auto& cand: cands.GetCandidates()){
		if (cand.second.GetIs3D()==1)++n3DBlobs;
		if (cand.second.GetClassifier()==goodBlob)++nGoodBlobs;

		//cout << "GOOD" << endl;	      
		int PID = cand.second.GetMCPID();
		int TopPID = cand.second.GetTopMCPID();
		int PTrackID = cand.second.GetMCParentTrackID();

		TVector3 FP = cand.second.GetFlightPath();
		double length = cand.second.GetDirection().Mag();
		double blobE = cand.second.GetTotalE();
		double dEdx = -1.0;
		if (length > 0.0) dEdx = blobE/length;
		double vtxDist = FP.Mag();
		double vtxZDist = abs(FP.Z());

		blobESum += blobE;
		/*
		if (PTrackID > nFSPart){
		  //Add something like this to learn how often this happened? ++nMultiIntBlobs;
		  continue;
		  }*/
		double candZ = cand.second.GetBegPos().Z();
		//if (cand.second.GetIs3D()==1) cout << "BlobIs3D" << endl;

		if (PTrackID==0 && !isPC){
		  map_hw_ALL_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_ALL_length_Recoil[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_blobE_Recoil[intType].univHist(universe)->Fill(blobE);
		  map_hw_ALL_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		  map_hw_ALL_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_ALL_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		  map_hw_ALL_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);

		  if (candZ > targetBoundary){
		    map_hw_tracker_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_tracker_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_blobE_Recoil[intType].univHist(universe)->Fill(blobE);
		    map_hw_tracker_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_tracker_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_tracker_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		    map_hw_tracker_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }

		  else{
		    map_hw_target_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_target_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_blobE_Recoil[intType].univHist(universe)->Fill(blobE);
		    map_hw_target_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_target_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_target_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		    map_hw_target_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }
		}

		else{
		  map_hw_ALL_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_ALL_length_Recoil[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_blobE_Recoil[intType].univHist(universe)->Fill(blobE);
		  map_hw_ALL_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		  map_hw_ALL_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_ALL_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		  map_hw_ALL_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);

		  if (candZ > targetBoundary){
		    map_hw_tracker_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_tracker_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_blobE_Recoil[intType].univHist(universe)->Fill(blobE);
		    map_hw_tracker_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_tracker_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_tracker_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		    map_hw_tracker_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }

		  else{
		    map_hw_target_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_target_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_blobE_Recoil[intType].univHist(universe)->Fill(blobE);
		    map_hw_target_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_target_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_target_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		    map_hw_target_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }
		}
	      }
	    }

	    //Event level Plots for passing Tejin Recoil but not Blob
	    if (leadBlobPasses) map_hw_leadBlob_passes_classifier_Recoil[intType].univHist(universe)->Fill(1);
	    else map_hw_leadBlob_passes_classifier_Recoil[intType].univHist(universe)->Fill(0);	  
	    map_hw_n3DBlobs_Recoil[intType].univHist(universe)->Fill(n3DBlobs);
	    map_hw_nGoodBlobs_Recoil[intType].univHist(universe)->Fill(nGoodBlobs);
	    map_hw_nBlobs_Recoil[intType].univHist(universe)->Fill(nBlobs);
	    map_hw_AvgBlobEnergy_Recoil[intType].univHist(universe)->Fill(blobESum/((double)(nBlobs)));
	    map_hw_RecoilEnergyGeV_Recoil[intType].univHist(universe)->Fill(recoilEnergy);
	  }
	  //Fails Tejin Recoil. I'm not going to treat the Tejin Blob Cut as special/independent of this recoil cut.
	  else {
	    NeutronCandidates::NeutCands cands = universe->GetCurrentNeutCands();
	    for (auto& cand: cands.GetCandidates()){
	      if (cand.second.GetIs3D()==1)++n3DBlobs;
	      if (cand.second.GetClassifier()==goodBlob)++nGoodBlobs;

	      //cout << "GOOD" << endl;	      
	      int PID = cand.second.GetMCPID();
	      int TopPID = cand.second.GetTopMCPID();
	      int PTrackID = cand.second.GetMCParentTrackID();

	      TVector3 FP = cand.second.GetFlightPath();
	      double length = cand.second.GetDirection().Mag();
	      double blobE = cand.second.GetTotalE();
	      double dEdx = -1.0;
	      if (length > 0.0) dEdx = blobE/length;
	      double vtxDist = FP.Mag();
	      double vtxZDist = abs(FP.Z());

	      blobESum += blobE;
	      /*
	      if (PTrackID > nFSPart){
		//Add something like this to learn how often this happened? ++nMultiIntBlobs;
		continue;
		}*/
	      double candZ = cand.second.GetBegPos().Z();
	      //if (cand.second.GetIs3D()==1) cout << "BlobIs3D" << endl;

	      if (PTrackID==0 && !isPC){
		map_hw_ALL_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		map_hw_ALL_length_CCQE[intType].univHist(universe)->Fill(length);
		map_hw_ALL_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		map_hw_ALL_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		map_hw_ALL_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		map_hw_ALL_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);

		if (candZ > targetBoundary){
		  map_hw_tracker_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_tracker_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_tracker_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_tracker_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		  map_hw_tracker_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_tracker_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		}

		else{ 
		  map_hw_target_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_target_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_target_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_target_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		  map_hw_target_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_target_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		}
	      }

	      else{
		map_hw_ALL_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		map_hw_ALL_length_CCQE[intType].univHist(universe)->Fill(length);
		map_hw_ALL_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		map_hw_ALL_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		map_hw_ALL_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		map_hw_ALL_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);

		if (candZ > targetBoundary){
		  map_hw_tracker_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_tracker_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_tracker_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_tracker_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		  map_hw_tracker_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_tracker_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		}

		else{
		  map_hw_target_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_target_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_target_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_target_blobE_CCQE[intType].univHist(universe)->Fill(blobE);
		  map_hw_target_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_target_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		}
	      }
	    }
	  }

	  //Event Level Plots for passing Only the CCQE level no Recoil
	  if (leadBlobPasses) map_hw_leadBlob_passes_classifier_CCQE[intType].univHist(universe)->Fill(1);
	  else map_hw_leadBlob_passes_classifier_CCQE[intType].univHist(universe)->Fill(0);	  
	  map_hw_n3DBlobs_CCQE[intType].univHist(universe)->Fill(n3DBlobs);
	  map_hw_nGoodBlobs_CCQE[intType].univHist(universe)->Fill(nGoodBlobs);
	  map_hw_nBlobs_CCQE[intType].univHist(universe)->Fill(nBlobs);
	  map_hw_AvgBlobEnergy_CCQE[intType].univHist(universe)->Fill(blobESum/((double)(nBlobs)));
	  map_hw_RecoilEnergyGeV_CCQE[intType].univHist(universe)->Fill(recoilEnergy);
	}
      }
    }
  }

  TFile* outFile = new TFile((TString)(outDir)+"runEventLoop_sample_"+sampleNames[sample]+"_region_"+regionNames[region]+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+TString(to_string(nEntries))+"_Events.root","RECREATE");
  cout << "Writing" << endl;
  for (auto band : error_bands){
    int i=0;
    vector<CVUniverse*> error_band_universes = band.second;
    for (auto universe : error_band_universes){
      Write1DHistsToFile(histsALL, universe, outFile);
    }
  }

  outFile->Close();
  cout << "HEY YOU DID IT!!!" << endl;
  return 0;

}
