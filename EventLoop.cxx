//File: EventLoop.cxx
//Info: This is a script to run a loop over all events in a single nTuple file and perform some plotting. Will eventually exist as the basis for the loops over events in analysis.
//
//Usage: EventLoop.cxx <single_MasterAnaDev_NTuple_input>
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

bool PassesTgt5FVCuts(CVUniverse& univ){
  vector<double> vtx = univ.GetVtx();
  double side = 850.0*2.0/sqrt(3.0);
  if (vtx[2] < 5756.71 || vtx[2] > 5801.24 ) return false;
  if (fabs(vtx[0]) > 850.00) return false;
  double slope = 1.0/sqrt(3.0);
  
  if(fabs(vtx[1]) < side - slope*fabs(vtx[0]))return true;
  return false;
}

bool PassesCleanCCAntiNuCuts(CVUniverse& univ){
  return 
    (univ.GetNDeadDiscriminatorsUpstreamMuon() < 2) &&
    (univ.GetNuHelicity() == 2) &&
    //    (univ.GetIsMinosMatchTrackOLD() == 1) &&
    (univ.GetIsMinosMatchTrack() == 1) &&
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

bool PassesCuts(CVUniverse& univ){
  return
    PassesTgt5FVCuts(univ) &&
    PassesCleanCCAntiNuCuts(univ) &&
    PassesTejinCCQECuts(univ);
}

int main(int argc, char* argv[]) {

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc != 2) {
    cout << "You don't understand this do you... You need a single input file!!!!!!!!!!!" << endl;
    return 1;
  }
  
  PlotUtils::ChainWrapper* chain = makeChainWrapperPtr(string(argv[1]),"MasterAnaDev");
  
  CVUniverse* CV = new CVUniverse(chain);
  map< string, vector<CVUniverse*>> error_bands;
  error_bands[string("CV")].push_back(CV);
  
  PlotUtils::HistWrapper<CVUniverse> hw_nBlobs("hw_nBlobs","Number of Neutron Blobs for this selection MAD 6D \"Full\" (04-28-2021 10:11 AM)",100,0.0,100.0,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_nBlobs_from_CandObj("hw_nBlobs_fromCandObj","Number of Neutron Blobs for this selection MAD 6D \"Full\" (04-28-2021 10:11 AM)",100,0.0,100.0,error_bands);
  
  int nEntries = chain->GetEntries();
  cout << "Processing " << nEntries << " events." << endl;
  for (int i=0; i<nEntries;++i){
    if (i%(nEntries/100)==0) cout << (100*i)/nEntries << "% finished." << endl;
    for (auto band : error_bands){
      vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes){
	universe->SetEntry(i);
	universe->UpdateNeutCands();
	if (PassesCuts(*universe)){
	  hw_nBlobs.univHist(universe)->Fill(universe->GetNNeutBlobs());
	  hw_nBlobs_from_CandObj.univHist(universe)->Fill(universe->GetNNeutCands());
	}
      }
    }
  }

  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->cd();
  for (auto band : error_bands){
    int i=0;
    vector<CVUniverse*> error_band_universes = band.second;
    for (auto universe : error_band_universes){
      ++i;
      hw_nBlobs.univHist(universe)->Draw();
      c1->Print("/minerva/app/users/dlast/TEST_MAT_Plots/h_nBlobs_REFERENCE_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_neutMC.pdf");
      hw_nBlobs_from_CandObj.univHist(universe)->Draw();
      c1->Print("/minerva/app/users/dlast/TEST_MAT_Plots/h_nBlobs_from_CandObj_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_neutMCtest.pdf");
    }
  }

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;

}
