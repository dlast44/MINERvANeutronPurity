//File: All1DIntTypeStackedPlots.cxx
//Info: This is a script to run a loop over all MC int type sorted plots in a single histos file and save nice plots from them.
//
//Usage: All1DIntTypeStackedPlots <histos_file> <output_directory> <name_of_data_sample>
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
#include "TKey.h"
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
#include "PlotUtils/MnvH1D.h"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

using namespace std;
using namespace PlotUtils;

TCanvas* DrawToCanvas(string name_QE, TFile* file, TString sample){

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->cd();

  MnvH1D* h_QE = (MnvH1D*)file->Get((TString)name_QE);
  h_QE->SetLineColor(kBlue);
  h_QE->SetFillColor(kBlue);

  string name = (string)h_QE->GetName();
  name.erase(name.length()-3,name.length());
  cout << "Handling: " << name << endl;
  string title = (string)h_QE->GetTitle();
  TString Xtitle = h_QE->GetXaxis()->GetTitle();
  TString Ytitle = h_QE->GetYaxis()->GetTitle();
  cout << title << endl;
  title.erase(0, 8);
  cout << title << endl;
  cout << "" << endl;

  MnvH1D* h_RES = (MnvH1D*)file->Get((TString)name+"_RES");
  h_RES->SetLineColor(kRed);
  h_RES->SetFillColor(kRed);

  MnvH1D* h_DIS = (MnvH1D*)file->Get((TString)name+"_DIS");
  h_DIS->SetLineColor(kGreen);
  h_DIS->SetFillColor(kGreen);

  MnvH1D* h_2p2h = (MnvH1D*)file->Get((TString)name+"_2p2h");
  h_2p2h->SetLineColor(kBlack);
  h_2p2h->SetFillColor(kBlack);

  MnvH1D* h_Other = (MnvH1D*)file->Get((TString)name+"_Other");
  h_Other->SetLineColor(kGray);
  h_Other->SetFillColor(kGray);

  THStack* h = new THStack();
  h->Add(h_Other);
  h->Add(h_2p2h);
  h->Add(h_DIS);
  h->Add(h_RES);
  h->Add(h_QE);

  h->Draw("hist");
  c1->Update();

  h->SetTitle(sample+" "+title.c_str());
  h->GetXaxis()->SetTitle(Xtitle);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitle(Ytitle);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleOffset(1.075);
  
  size_t pos=0;
  if ((pos=name.find("_primary_parent")) != string::npos){
    h->GetXaxis()->SetBinLabel(1,"Other");
    h->GetXaxis()->SetBinLabel(2,"");
    h->GetXaxis()->SetBinLabel(3,"n");
    h->GetXaxis()->SetBinLabel(4,"p");
    h->GetXaxis()->SetBinLabel(5,"#pi^{0}");
    h->GetXaxis()->SetBinLabel(6,"#pi^{+}");
    h->GetXaxis()->SetBinLabel(7,"#pi^{-}");
    h->GetXaxis()->SetBinLabel(8,"#gamma");
    h->GetXaxis()->SetBinLabel(9,"e");
    h->GetXaxis()->SetBinLabel(10,"#mu");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitle("Blob Primary Parent");
    h->GetXaxis()->SetTitleSize(0.045);
  }

  h->Draw("hist");
  c1->Update();

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  leg->AddEntry(h_QE,"QE");
  leg->AddEntry(h_RES,"RES");
  leg->AddEntry(h_DIS,"DIS");
  leg->AddEntry(h_2p2h,"2p2h");
  leg->AddEntry(h_Other,"Other");

  leg->Draw();
  c1->Update();
  return c1;
}

bool PathExists(string path){
  struct stat buffer;
  return (stat (path.c_str(), &buffer) == 0);
}

int main(int argc, char* argv[]) {

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc != 4) {
    cout << "Check usage..." << endl;
    return 2;
  }

  string inName=string(argv[1]);
  string outDir=string(argv[2]);
  TString sample = argv[3];

  if (PathExists(outDir)){
    cout << "Thank you for choosing a path for output files that exists." << endl;
  }
  else{
    cout << "Output directory doesn't exist. Exiting" << endl;
    return 3;
  }

  string rootExt = ".root";
  string slash = "/";
  string token;
  string inNameStub = inName;
  size_t pos=0;

  //cout << inNameStub << endl;
  while ((pos = inNameStub.find(slash)) != string::npos){
    //cout << inNameStub << endl;
    token = inNameStub.substr(0, pos);
    //cout << token << endl;
    inNameStub.erase(0, pos+slash.length());
  }
  //cout << inNameStub << endl;
  if ((pos=inNameStub.find(rootExt)) == string::npos){
    cout << "Input need be .root file." << endl;
    return 4;
  }

  cout << "Input file name parsed to: " << inNameStub << endl;

  TFile* inFile = new TFile(inName.c_str(),"READ");

  TList* keyList = inFile->GetListOfKeys();
  if (!keyList){
    cout << "List of keys failed to get." << endl;
    return 5;
  }
  TIter next(keyList);
  TKey* key;
  while ( key = (TKey*)next() ){
    //cout << key->GetName() << endl;
    pos=0;
    string name=(string)key->GetName();
    if((pos = name.find("_QE")) == string::npos) continue;
    TCanvas* c1 = DrawToCanvas(name,inFile,sample);
    name.erase(name.length()-3,name.length());
    c1->Print((TString)outDir+(TString)name+"_stacked_test.pdf");
    c1->Print((TString)outDir+(TString)name+"_stacked_test.png");
    delete c1;
  }

  //parsingTest("name",inFile);

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
