//File: EventLoop.cxx
//Info: This is a script to run a loop over all events in a single nTuple file and perform some plotting. Will eventually exist as the basis for the loops over events in analysis.
//
//Usage: EventLoop.cxx <MasterAnaDev_NTuple_list/single_file> <0=MC/1=PC> <0=tracker/1=targets/2=both> <output_directory> <tag_for_naming_files> optional: <n_event>
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

bool PassesCleanCCAntiNuCuts(CVUniverse& univ, int isPC){
  int MINOSMatch=0;
  if (isPC){
    //    MINOSMatch=univ.GetIsMinosMatchTrackOLD();
    MINOSMatch=1;
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
bool PassesTejinRecoilCut(CVUniverse& univ, int isPC){
  if (isPC) return true;
  double Q2GeV = univ.GetQ2QEPickledGeV();
  double recoilEGeV = univ.GetRecoilEnergyGeV();
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
    (univ.GetNImprovedMichel() > 0);
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

bool PassesCuts(CVUniverse& univ, int isPC, int region){
  //cout << "HELLO" << endl;
  return
    PassesFVCuts(univ, region) &&
    PassesCleanCCAntiNuCuts(univ, isPC) &&
    PassesTejinCCQECuts(univ);
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
  if (argc < 6 || argc > 7) {
    cout << "Check usage..." << endl;
    return 2;
  }

  string playlist=string(argv[1]);
  int isPC=atoi(argv[2]);
  int region=atoi(argv[3]);
  string outDir=string(argv[4]);
  string tag=string(argv[5]);
  int nEntries=0;

  if (argc == 7){
    nEntries=atoi(argv[6]);
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

  for(auto type: typeNames){
    map_hw_tracker_primary_parent_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_primary_parent_CCQE_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE)",10,0,10,error_bands);
    map_hw_tracker_primary_parent_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_primary_parent_Recoil_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil)",10,0,10,error_bands);
    map_hw_tracker_primary_parent_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_primary_parent_Tejin_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil, Blob)",10,0,10,error_bands);
    map_hw_tracker_primary_parent_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_primary_parent_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil, Blob, Tracker Blob)",10,0,10,error_bands);

    map_hw_tracker_length_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_length_CCQE_"+type.second,"True "+type.second+" Blob Length (CCQE)",50,0,500,error_bands);
    map_hw_tracker_length_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_length_Recoil_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil)",50,0,500,error_bands);
    map_hw_tracker_length_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_length_Tejin_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil, Blob)",50,0,500,error_bands);
    map_hw_tracker_length_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_length_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil, Blob, Tracker Blob)",50,0,500,error_bands);

    map_hw_tracker_avg_dEdx_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_avg_dEdx_CCQE_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE)",50,0,50,error_bands);
    map_hw_tracker_avg_dEdx_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_avg_dEdx_Recoil_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil)",50,0,50,error_bands);
    map_hw_tracker_avg_dEdx_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_avg_dEdx_Tejin_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil, Blob)",50,0,50,error_bands);
    map_hw_tracker_avg_dEdx_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_avg_dEdx_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil, Blob, Tracker Blob)",50,0,50,error_bands);

    map_hw_tracker_dist_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_dist_CCQE_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE)",300,0,3000,error_bands);
    map_hw_tracker_dist_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_dist_Recoil_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil)",300,0,3000,error_bands);
    map_hw_tracker_dist_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_dist_Tejin_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil, Blob)",300,0,3000,error_bands);
    map_hw_tracker_dist_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_dist_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil, Blob, Tracker Blob)",300,0,3000,error_bands);

    map_hw_tracker_Zdist_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_Zdist_CCQE_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE)",300,0,3000,error_bands);
    map_hw_tracker_Zdist_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_Zdist_Recoil_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil)",300,0,3000,error_bands);
    map_hw_tracker_Zdist_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_Zdist_Tejin_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil, Blob)",300,0,3000,error_bands);
    map_hw_tracker_Zdist_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_tracker_Zdist_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil, Blob, Tracker Blob)",300,0,3000,error_bands);

    map_hw_target_primary_parent_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_primary_parent_CCQE_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE)",10,0,10,error_bands);
    map_hw_target_primary_parent_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_primary_parent_Recoil_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil)",10,0,10,error_bands);
    map_hw_target_primary_parent_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_primary_parent_Tejin_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil, Blob)",10,0,10,error_bands);
    map_hw_target_primary_parent_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_primary_parent_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil, Blob, Tracker Blob)",10,0,10,error_bands);

    map_hw_target_length_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_length_CCQE_"+type.second,"True "+type.second+" Blob Length (CCQE)",50,0,500,error_bands);
    map_hw_target_length_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_length_Recoil_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil)",50,0,500,error_bands);
    map_hw_target_length_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_length_Tejin_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil, Blob)",50,0,500,error_bands);
    map_hw_target_length_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_length_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil, Blob, Tracker Blob)",50,0,500,error_bands);

    map_hw_target_avg_dEdx_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_avg_dEdx_CCQE_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE)",50,0,50,error_bands);
    map_hw_target_avg_dEdx_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_avg_dEdx_Recoil_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil)",50,0,50,error_bands);
    map_hw_target_avg_dEdx_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_avg_dEdx_Tejin_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil, Blob)",50,0,50,error_bands);
    map_hw_target_avg_dEdx_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_avg_dEdx_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil, Blob, Tracker Blob)",50,0,50,error_bands);

    map_hw_target_dist_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_dist_CCQE_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE)",300,0,3000,error_bands);
    map_hw_target_dist_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_dist_Recoil_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil)",300,0,3000,error_bands);
    map_hw_target_dist_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_dist_Tejin_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil, Blob)",300,0,3000,error_bands);
    map_hw_target_dist_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_dist_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil, Blob, Tracker Blob)",300,0,3000,error_bands);

    map_hw_target_Zdist_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_Zdist_CCQE_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE)",300,0,3000,error_bands);
    map_hw_target_Zdist_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_Zdist_Recoil_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil)",300,0,3000,error_bands);
    map_hw_target_Zdist_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_Zdist_Tejin_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil, Blob)",300,0,3000,error_bands);
    map_hw_target_Zdist_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_target_Zdist_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil, Blob, Tracker Blob)",300,0,3000,error_bands);

    map_hw_ALL_primary_parent_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_primary_parent_CCQE_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE)",10,0,10,error_bands);
    map_hw_ALL_primary_parent_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_primary_parent_Recoil_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil)",10,0,10,error_bands);
    map_hw_ALL_primary_parent_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_primary_parent_Tejin_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil, Blob)",10,0,10,error_bands);
    map_hw_ALL_primary_parent_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_primary_parent_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Primary Particle Matched To Blob (CCQE, Recoil, Blob, Tracker Blob)",10,0,10,error_bands);

    map_hw_ALL_length_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_length_CCQE_"+type.second,"True "+type.second+" Blob Length (CCQE)",50,0,500,error_bands);
    map_hw_ALL_length_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_length_Recoil_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil)",50,0,500,error_bands);
    map_hw_ALL_length_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_length_Tejin_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil, Blob)",50,0,500,error_bands);
    map_hw_ALL_length_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_length_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Length (CCQE, Recoil, Blob, Tracker Blob)",50,0,500,error_bands);

    map_hw_ALL_avg_dEdx_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_avg_dEdx_CCQE_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE)",50,0,50,error_bands);
    map_hw_ALL_avg_dEdx_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_avg_dEdx_Recoil_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil)",50,0,50,error_bands);
    map_hw_ALL_avg_dEdx_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_avg_dEdx_Tejin_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil, Blob)",50,0,50,error_bands);
    map_hw_ALL_avg_dEdx_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_avg_dEdx_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Energy/Length (CCQE, Recoil, Blob, Tracker Blob)",50,0,50,error_bands);

    map_hw_ALL_dist_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_dist_CCQE_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE)",300,0,3000,error_bands);
    map_hw_ALL_dist_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_dist_Recoil_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil)",300,0,3000,error_bands);
    map_hw_ALL_dist_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_dist_Tejin_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil, Blob)",300,0,3000,error_bands);
    map_hw_ALL_dist_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_dist_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Dist. To Vtx. (CCQE, Recoil, Blob, Tracker Blob)",300,0,3000,error_bands);

    map_hw_ALL_Zdist_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_Zdist_CCQE_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE)",300,0,3000,error_bands);
    map_hw_ALL_Zdist_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_Zdist_Recoil_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil)",300,0,3000,error_bands);
    map_hw_ALL_Zdist_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_Zdist_Tejin_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil, Blob)",300,0,3000,error_bands);
    map_hw_ALL_Zdist_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_ALL_Zdist_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Blob Absolute Z Dist. To Vtx. (CCQE, Recoil, Blob, Tracker Blob)",300,0,3000,error_bands);

    map_hw_leadBlob_passes_classifier_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_leadBlob_passes_classifier_CCQE_"+type.second,"True "+type.second+" Leading Blob Passes (CCQE)",2,0,2,error_bands);
    map_hw_leadBlob_passes_classifier_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_leadBlob_passes_classifier_Recoil_"+type.second,"True "+type.second+" Leading Blob Passes (CCQE, Recoil)",2,0,2,error_bands);
    map_hw_leadBlob_passes_classifier_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_leadBlob_passes_classifier_Tejin_"+type.second,"True "+type.second+" Leading Blob Passes (CCQE, Recoil, Blob)",2,0,2,error_bands);
    map_hw_leadBlob_passes_classifier_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_leadBlob_passes_classifier_Tejin_TrackerONLY_"+type.second,"True "+type.second+" Leading Blob Passes (CCQE, Recoil, Blob, Tracker Blob)",2,0,2,error_bands);

    map_hw_n3DBlobs_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_n3DBlobs_CCQE_"+type.second,"True "+type.second+" No. 3D Blobs (CCQE)",10,0,10,error_bands);
    map_hw_n3DBlobs_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_n3DBlobs_Recoil_"+type.second,"True "+type.second+" No. 3D Blobs (CCQE, Recoil)",10,0,10,error_bands);
    map_hw_n3DBlobs_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_n3DBlobs_Tejin_"+type.second,"True "+type.second+" No. 3D Blobs (CCQE, Recoil, Blob)",10,0,10,error_bands);
    map_hw_n3DBlobs_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_n3DBlobs_Tejin_TrackerONLY_"+type.second,"True "+type.second+" No. 3D Blobs (CCQE, Recoil, Blob, Tracker Blob)",10,0,10,error_bands);

    map_hw_nGoodBlobs_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nGoodBlobs_CCQE_"+type.second,"True "+type.second+" No. Blobs Which Pass (CCQE)",10,0,10,error_bands);
    map_hw_nGoodBlobs_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nGoodBlobs_Recoil_"+type.second,"True "+type.second+" No. Blobs Which Pass (CCQE, Recoil)",10,0,10,error_bands);
    map_hw_nGoodBlobs_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nGoodBlobs_Tejin_"+type.second,"True "+type.second+" No. Blobs Which Pass (CCQE, Recoil, Blob)",10,0,10,error_bands);
    map_hw_nGoodBlobs_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nGoodBlobs_Tejin_TrackerONLY_"+type.second,"True "+type.second+" No. Blobs Which Pass (CCQE, Recoil, Blob, Tracker Blob)",10,0,10,error_bands);

    map_hw_nBlobs_CCQE[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nBlobs_CCQE_"+type.second,"True "+type.second+" No. of Blobs (CCQE)",100,0,100,error_bands);
    map_hw_nBlobs_Recoil[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nBlobs_Recoil_"+type.second,"True "+type.second+" No. of Blobs (CCQE, Recoil)",100,0,100,error_bands);
    map_hw_nBlobs_Tejin[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nBlobs_Tejin_"+type.second,"True "+type.second+" No. of Blobs (CCQE, Recoil, Blob)",100,0,100,error_bands);
    map_hw_nBlobs_Tejin_TrackerONLY[type.first]=PlotUtils::HistWrapper<CVUniverse>("hw_nBlobs_Tejin_TrackerONLY_"+type.second,"True "+type.second+" No. of Blobs (CCQE, Recoil, Blob, Tracker Blob)",100,0,100,error_bands);
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
    &map_hw_nBlobs_Tejin_TrackerONLY[1],&map_hw_nBlobs_Tejin_TrackerONLY[2],&map_hw_nBlobs_Tejin_TrackerONLY[3],&map_hw_nBlobs_Tejin_TrackerONLY[8],&map_hw_nBlobs_Tejin_TrackerONLY[0]
  };

  if(!nEntries) nEntries = chain->GetEntries();
  cout << "Processing " << nEntries << " events." << endl;
  int n3DBlobs=0;
  int nGoodBlobs=0;
  for (int i=0; i<nEntries;++i){
    if (i%(nEntries/100)==0) cout << (100*i)/nEntries << "% finished." << endl;
    //if (i%(10000)==0) cout << i << " entries finished." << endl;
    for (auto band : error_bands){
      vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes){
	n3DBlobs=0;
	nGoodBlobs=0;
	universe->SetEntry(i);
	universe->UpdateNeutCands();
	int nBlobs = universe->GetNNeutCands();
	
	//Passes CCQE Cuts that matche Tejin's selection
	if (PassesCuts(*universe, isPC, region)){
	  
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
	  if (PassesTejinRecoilCut(*universe, isPC)){
	    
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
		if (length > 0.0) dEdx = cand.second.GetTotalE()/length;
		double vtxDist = FP.Mag();
		double vtxZDist = abs(FP.Z());

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
		    map_hw_ALL_dist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_ALL_Zdist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxZDist);

		    if (candZ > targetBoundary){
		      map_hw_tracker_primary_parent_Tejin_TrackerONLY[intType].univHist(universe)->Fill(PDGbins[PID]);
		      map_hw_tracker_length_Tejin_TrackerONLY[intType].univHist(universe)->Fill(length);
		      map_hw_tracker_avg_dEdx_Tejin_TrackerONLY[intType].univHist(universe)->Fill(dEdx);
		      map_hw_tracker_dist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxDist);
		      map_hw_tracker_Zdist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxZDist);
		    }

		    else {
		      map_hw_target_primary_parent_Tejin_TrackerONLY[intType].univHist(universe)->Fill(PDGbins[PID]);
		      map_hw_target_length_Tejin_TrackerONLY[intType].univHist(universe)->Fill(length);
		      map_hw_target_avg_dEdx_Tejin_TrackerONLY[intType].univHist(universe)->Fill(dEdx);
		      map_hw_target_dist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxDist);
		      map_hw_target_Zdist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxZDist);
		    }
		  }

		  map_hw_ALL_primary_parent_Tejin[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_ALL_length_Tejin[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_Tejin[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_dist_Tejin[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_Tejin[intType].univHist(universe)->Fill(vtxZDist);
		  map_hw_ALL_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_ALL_length_Recoil[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		  map_hw_ALL_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_ALL_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);

		  if (candZ > targetBoundary){
		    map_hw_tracker_primary_parent_Tejin[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_tracker_length_Tejin[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_Tejin[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_dist_Tejin[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_Tejin[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_tracker_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_tracker_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_tracker_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_tracker_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }

		  else{
		    map_hw_target_primary_parent_Tejin[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_target_length_Tejin[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_Tejin[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_dist_Tejin[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_Tejin[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_target_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_target_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_target_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_target_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
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
		    map_hw_ALL_dist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_ALL_Zdist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxZDist);

		    if (candZ > targetBoundary){
		      map_hw_tracker_primary_parent_Tejin_TrackerONLY[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		      map_hw_tracker_length_Tejin_TrackerONLY[intType].univHist(universe)->Fill(length);
		      map_hw_tracker_avg_dEdx_Tejin_TrackerONLY[intType].univHist(universe)->Fill(dEdx);
		      map_hw_tracker_dist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxDist);
		      map_hw_tracker_Zdist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxZDist);
		    }		  

		    else {
		      map_hw_target_primary_parent_Tejin_TrackerONLY[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		      map_hw_target_length_Tejin_TrackerONLY[intType].univHist(universe)->Fill(length);
		      map_hw_target_avg_dEdx_Tejin_TrackerONLY[intType].univHist(universe)->Fill(dEdx);
		      map_hw_target_dist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxDist);
		      map_hw_target_Zdist_Tejin_TrackerONLY[intType].univHist(universe)->Fill(vtxZDist);
		    }
		  }

		  map_hw_ALL_primary_parent_Tejin[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_ALL_length_Tejin[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_Tejin[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_dist_Tejin[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_Tejin[intType].univHist(universe)->Fill(vtxZDist);
		  map_hw_ALL_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_ALL_length_Recoil[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		  map_hw_ALL_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_ALL_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);

		  if (candZ > targetBoundary){
		    map_hw_tracker_primary_parent_Tejin[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_tracker_length_Tejin[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_Tejin[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_dist_Tejin[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_Tejin[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_tracker_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_tracker_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_tracker_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_tracker_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }

		  else{
		    map_hw_target_primary_parent_Tejin[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_target_length_Tejin[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_Tejin[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_dist_Tejin[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_Tejin[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_target_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_target_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_target_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_target_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
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

	      if (TejinBlobValue==2){
		if (leadBlobPasses) map_hw_leadBlob_passes_classifier_Tejin_TrackerONLY[intType].univHist(universe)->Fill(1);
		else map_hw_leadBlob_passes_classifier_Tejin_TrackerONLY[intType].univHist(universe)->Fill(0);
		map_hw_n3DBlobs_Tejin_TrackerONLY[intType].univHist(universe)->Fill(n3DBlobs);
		map_hw_nGoodBlobs_Tejin_TrackerONLY[intType].univHist(universe)->Fill(nGoodBlobs);
		map_hw_nBlobs_Tejin_TrackerONLY[intType].univHist(universe)->Fill(nBlobs);
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
		double dEdx = -1.0;
		if (length > 0.0) dEdx = cand.second.GetTotalE()/length;
		double vtxDist = FP.Mag();
		double vtxZDist = abs(FP.Z());
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
		  map_hw_ALL_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		  map_hw_ALL_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_ALL_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);

		  if (candZ > targetBoundary){
		    map_hw_tracker_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_tracker_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_tracker_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_tracker_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }

		  else{
		    map_hw_target_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_target_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_target_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		    map_hw_target_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }
		}

		else{
		  map_hw_ALL_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_ALL_length_Recoil[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		  map_hw_ALL_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_ALL_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_ALL_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_ALL_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_ALL_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);

		  if (candZ > targetBoundary){
		    map_hw_tracker_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_tracker_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_tracker_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_tracker_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_tracker_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		    map_hw_tracker_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_tracker_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		  }

		  else{
		    map_hw_target_primary_parent_Recoil[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_target_length_Recoil[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_Recoil[intType].univHist(universe)->Fill(dEdx);
		    map_hw_target_dist_Recoil[intType].univHist(universe)->Fill(vtxDist);
		    map_hw_target_Zdist_Recoil[intType].univHist(universe)->Fill(vtxZDist);
		    map_hw_target_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		    map_hw_target_length_CCQE[intType].univHist(universe)->Fill(length);
		    map_hw_target_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
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
	      double dEdx = -1.0;
	      if (length > 0.0) dEdx = cand.second.GetTotalE()/length;
	      double vtxDist = FP.Mag();
	      double vtxZDist = abs(FP.Z());
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
		map_hw_ALL_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		map_hw_ALL_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);

		if (candZ > targetBoundary){
		  map_hw_tracker_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_tracker_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_tracker_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_tracker_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_tracker_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		}

		else{ 
		  map_hw_target_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[PID]);
		  map_hw_target_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_target_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_target_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_target_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		}
	      }

	      else{
		map_hw_ALL_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		map_hw_ALL_length_CCQE[intType].univHist(universe)->Fill(length);
		map_hw_ALL_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		map_hw_ALL_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		map_hw_ALL_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);

		if (candZ > targetBoundary){
		  map_hw_tracker_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_tracker_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_tracker_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
		  map_hw_tracker_dist_CCQE[intType].univHist(universe)->Fill(vtxDist);
		  map_hw_tracker_Zdist_CCQE[intType].univHist(universe)->Fill(vtxZDist);
		}

		else{
		  map_hw_target_primary_parent_CCQE[intType].univHist(universe)->Fill(PDGbins[TopPID]);
		  map_hw_target_length_CCQE[intType].univHist(universe)->Fill(length);
		  map_hw_target_avg_dEdx_CCQE[intType].univHist(universe)->Fill(dEdx);
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
	    
	}
      }
    }
  }

  TFile* outFile = new TFile((TString)(outDir)+"runEventLoop_region_"+regionNames[region]+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+TString(to_string(nEntries))+"_Events.root","RECREATE");
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
