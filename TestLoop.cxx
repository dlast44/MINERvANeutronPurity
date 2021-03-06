//File: TestLoop.cxx
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

  if(!nEntries) nEntries = chain->GetEntries();
  cout << "Processing " << nEntries << " events." << endl;
  for (int i=0; i<nEntries;++i){
    if (i%(nEntries/100)==0) cout << (100*i)/nEntries << "% finished." << endl;
    //if (i%(10000)==0) cout << i << " entries finished." << endl;
    for (auto band : error_bands){
      vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes){
	cout << "Event: " << i << endl;
	cout << "" << endl;
	universe->SetEntry(i);
	universe->UpdateNeutCands();
	for (auto& cand: universe->GetCurrentNeutCands().GetCandidates()){
	  cout << "ID: " << cand.second.GetID() << endl;
	  cout << "Z Dist: " << cand.second.GetFlightPath().Z() << endl;
	  cout << "TotalE: " << cand.second.GetTotalE() << endl;
	  cout << "Angle: " << cand.second.GetAngleToFP() << endl;
	  cout << "Is3D: " << cand.second.GetIs3D() << endl;
	  cout << "Classifier: " << cand.second.GetClassifier() << endl;
	  cout << "Direction Magnitude: " << cand.second.GetDirection().Mag() << endl;
	  cout << "Beg X: " << cand.second.GetBegPos().X() << endl;
	  cout << "End X: " << cand.second.GetEndPos().X() << endl;
	  cout << "Beg Y: " << cand.second.GetBegPos().Y() << endl;
	  cout << "End Y: " << cand.second.GetEndPos().Y() << endl;
	  cout << "Beg Z: " << cand.second.GetBegPos().Z() << endl;
	  cout << "End Z: " << cand.second.GetEndPos().Z() << endl;
	  cout << "" << endl;
	}
      }
    }	
  }

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;

}
