//File: testExecutable.cxx
//Info: This is a test executable for the CMake build structure to include ROOT and PlotUtils (and UnfoldUtils?) properly...
//
//Usage: testExecutable
//
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

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

int main(int argc, char* argv[]) {

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  if (argc > 1){
    std::cout << "You don't understand this do you..." << std::endl;
  }

  srand( (unsigned)time(NULL) );

  PlotUtils::MnvH1D* HelloWorld = new PlotUtils::MnvH1D("h_HelloWorld","Testing to see if I can properly handle a MnvH1D at the simplest level...",100,0.0,1.0);

  for (int i=0; i<100000;++i){
    HelloWorld->Fill((float)rand()/RAND_MAX);
  }

  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->cd();
  HelloWorld->Draw();
  c1->Print("HELLO_WORLD_WITH_CINTEX.pdf");

  std::cout << "HEY YOU DID IT" << std::endl;
  return 0;

}
