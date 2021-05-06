//File: EventLoop.cxx
//Info: This is a script to run a loop over all events in a single nTuple file and perform some plotting. Will eventually exist as the basis for the loops over events in analysis.
//
//Usage: EventLoop.cxx <MasterAnaDev_NTuple_list/single_file> <0=MC/1=PC>
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

bool PassesFVCuts(CVUniverse& univ){
  vector<double> vtx = univ.GetVtx();
  double side = 850.0*2.0/sqrt(3.0);
  if (vtx[2] > 5850.0 || vtx[2] < 4500.0 ) return false;
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
    MINOSMatch=univ.GetIsMinosMatchTrack();
  }
  return
    (univ.GetNDeadDiscriminatorsUpstreamMuon() < 2) &&
    (univ.GetNuHelicity() == 2) &&
    (MINOSMatch == 1) &&
    (TMath::RadToDeg()*univ.GetThetamu() < 20.0) &&
    (univ.GetPmu() < 20000.0 && univ.GetPmu() > 1500.0);
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

bool PassesTejinBlobCuts(CVUniverse& universe){
  NeutronCandidates::NeutCand leadingBlob=universe.GetCurrentLeadingNeutCand();
  int leading3D=leadingBlob.GetIs3D();
  TVector3 leadingFP=leadingBlob.GetFlightPath();
  TVector3 muonMom(universe.GetMuon4V().Z(),universe.GetMuon4V().Y(),universe.GetMuon4V().Z());
  if (leadingFP.Mag()==0 || muonMom.Mag()==0) return false;
  else{
    return
      (leading3D==1) &&
      (leadingFP.Angle(muonMom) > 0.261799388);
  }
}

bool PassesCuts(CVUniverse& univ, int isPC){
  //cout << "HELLO" << endl;
  return
    PassesFVCuts(univ) &&
    PassesCleanCCAntiNuCuts(univ, isPC) &&
    PassesTejinCCQECuts(univ);
}

int main(int argc, char* argv[]) {

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc != 3) {
    cout << "You don't understand this do you... You need a single input file and !!!!!!!!!!!" << endl;
    return 1;
  }

  int isPC=atoi(argv[2]);
  
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

  PlotUtils::ChainWrapper* chain = makeChainWrapperPtr(string(argv[1]),"MasterAnaDev");
  
  CVUniverse* CV = new CVUniverse(chain);
  map<string, vector<CVUniverse*>> error_bands;
  error_bands[string("CV")].push_back(CV);

  int nEntries = chain->GetEntries();
  cout << "Processing " << nEntries << " events." << endl;
  for (int i=0; i<nEntries;++i){
    if (i%(nEntries/100)==0) cout << (100*i)/nEntries << "% finished." << endl;
    for (auto band : error_bands){
      vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes){
	universe->SetEntry(i);
	universe->UpdateNeutCands();
	if (PassesCuts(*universe, isPC)){
	  NeutronCandidates::NeutCand leadingBlob=universe->GetCurrentLeadingNeutCand();
	  int is3D = leadingBlob.GetIs3D();
	  double ang = leadingBlob.GetAngleToFP();
	  TVector3 FP = leadingBlob.GetFlightPath();
	  TVector3 dir = leadingBlob.GetDirection();
	  TVector3 pos = leadingBlob.GetBegPos();
	  TVector3 end = leadingBlob.GetEndPos();
	  TVector3 vtx = leadingBlob.GetEvtVtx();
	  if (is3D==1){
	    cout << "" << endl;
	    cout << "Beg Position: " << pos.X() << " " << pos.Y() << " " << pos.Z() << endl;
	    cout << "End Position: " << end.X() << " " << end.Y() << " " << end.Z() << endl;
	    cout << "Blob Direction: " << dir.X() << " " << dir.Y() << " " << dir.Z() << endl;
	    cout << "Vtx: " << vtx.X() << " " << vtx.Y() << " " << vtx.Z() << endl;
	    cout << "Flight Path: " << FP.X() << " " << FP.Y() << " " << FP.Z() << endl;
	    cout << "Angle: " << ang << endl;
	  }
	}
      }
    }
  }

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;

}
