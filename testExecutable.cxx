//File: testExecutable.cxx
//Info: This is a test executable for the CMake build structure to include ROOT and PlotUtils (and UnfoldUtils?) properly...
//
//Usage: testExecutable <single_MasterAnaDev_NTuple_input>
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

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

bool PassesCuts(CVUniverse& univ){
  bool tmp = univ.SetNFluxUniverses(0);
  return tmp;
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
  std::map< std::string, std::vector<CVUniverse*>> error_bands;
  error_bands["CV"].push_back(CV);
  
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
  hw_nTracks.hist->Draw();
  c1->Print("/minerva/app/users/dlast/BASIC_MAT_nTracks_Demonstrations3.pdf");
  hw_nBlobs.hist->Draw();
  c1->Print("/minerva/app/users/dlast/BASIC_MAT_nBlobs_Demonstrations3.pdf");

  std::cout << "HEY YOU DID IT!!!" << std::endl;
  return 0;

}
