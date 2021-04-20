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

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"

#include "syst/CVUniverse.h"
//#include "syst/NTracksShiftUniverse.h"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

bool PassesCuts(CVUniverse& univ){
  return univ.GetNTracks() > 0;
}

int main(int argc, char* argv[]) {

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc != 2) {
    std::cout << "You don't understand this do you... You need a single input file!!!!!!!!!!!" << std::endl;
    return 1;
  }
  
  PlotUtils::ChainWrapper* chain = makeChainWrapperPtr(std::string(argv[1]),"MasterAnaDev");
  
  CVUniverse* CV = new CVUniverse(chain);
  //NTracksShiftUniverse* NT_P1 = new NTracksShiftUniverse(chain, +1);
  //NTracksShiftUniverse* NT_M1 = new NTracksShiftUniverse(chain, -1);
  std::map< std::string, std::vector<CVUniverse*>> error_bands;
  error_bands[std::string("CV")].push_back(CV);
  //error_bands[std::string("NT")].push_back(NT_P1);
  //error_bands[std::string("NT")].push_back(NT_M1);
  
  PlotUtils::HistWrapper<CVUniverse> hw_nTracks("hw_nTracks","Check this against the input file "+TString(argv[1])+" multiplicity",10,0.0,10.0,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_nBlobs("hw_nBlobs","Check this against the input file "+TString(argv[1])+" MasterAnaDev_BlobIs3D_sz",100,0.0,100.0,error_bands);
  //PlotUtils::MnvH1D* h_nBlobs = new PlotUtils::MnvH1D("h_nBlobs","Check this against the input file "+TString(argv[1])+" MasterAnaDev_BlobIs3D_sz",100,0.0,100.0);
  
  for (int i=0; i<chain->GetEntries();++i){
    for (auto band : error_bands){
      std::vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes){
	universe->SetEntry(i);
	if (PassesCuts(*universe)){
	  hw_nTracks.univHist(universe)->Fill(universe->GetNTracks());
	  hw_nBlobs.univHist(universe)->Fill(universe->GetNNeutBlobs());
	}
      }
    }
  }

  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->cd();
  for (auto band : error_bands){
    int i=0;
    std::vector<CVUniverse*> error_band_universes = band.second;
    for (auto universe : error_band_universes){
      ++i;
      hw_nTracks.univHist(universe)->Draw();
      c1->Print("/minerva/app/users/dlast/TEST_MAT_Plots/"+(TString)(universe->ShortName())+(TString)(std::to_string(i))+".pdf");
    }
  }

  std::cout << "HEY YOU DID IT!!!" << std::endl;
  return 0;

}
