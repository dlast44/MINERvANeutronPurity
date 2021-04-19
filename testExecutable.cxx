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
#include "PlotUtile/MnvH1D"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

//UnfoldUtils includes???
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include "MinervaUnfold/MnvUnfold.h"

int main(int argc, char* argv[]) {
  if (argc > 1){
    std::cout << "You don't understand this do you..." << std::endl;
  }

  std::cout << "HEY YOU DID IT" << std::endl;
  return 0;

}
