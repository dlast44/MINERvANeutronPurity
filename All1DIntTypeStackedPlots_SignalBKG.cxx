//File: All1DIntTypeStackedPlots_SignalBKG.cxx
//Info: This is a script to run a loop over all MC int type sorted plots in a single histos file and save nice plots from them.
//
//Usage: All1DIntTypeStackedPlots_SignalBKG <histos_file_signal> <histos_file_BKG> <output_directory> <data_sample_name>
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

TCanvas* DrawToCanvas(string name_QE, TFile* sig_file, TFile* bkg_file, TString sample){

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->cd();

  MnvH1D* h_QE_Sig = (MnvH1D*)sig_file->Get((TString)name_QE);
  h_QE_Sig->SetLineColor(kBlue);
  h_QE_Sig->SetFillColor(kBlue);

  string name = (string)h_QE_Sig->GetName();
  name.erase(name.length()-3,name.length());
  cout << "Handling: " << name << endl;
  string title = (string)h_QE_Sig->GetTitle();
  TString Xtitle = h_QE_Sig->GetXaxis()->GetTitle();
  TString Ytitle = h_QE_Sig->GetYaxis()->GetTitle();
  cout << title << endl;
  title.erase(0, 8);
  cout << title << endl;
  cout << "" << endl;

  MnvH1D* h_RES_Sig = (MnvH1D*)sig_file->Get((TString)name+"_RES");
  h_RES_Sig->SetLineColor(kRed);
  h_RES_Sig->SetFillColor(kRed);

  MnvH1D* h_DIS_Sig = (MnvH1D*)sig_file->Get((TString)name+"_DIS");
  h_DIS_Sig->SetLineColor(kGreen);
  h_DIS_Sig->SetFillColor(kGreen);

  MnvH1D* h_2p2h_Sig = (MnvH1D*)sig_file->Get((TString)name+"_2p2h");
  h_2p2h_Sig->SetLineColor(kBlack);
  h_2p2h_Sig->SetFillColor(kBlack);

  MnvH1D* h_Other_Sig = (MnvH1D*)sig_file->Get((TString)name+"_Other");
  h_Other_Sig->SetLineColor(kGray);
  h_Other_Sig->SetFillColor(kGray);

  MnvH1D* h_QE_Bkg = (MnvH1D*)bkg_file->Get((TString)name+"_QE");
  h_QE_Bkg->SetLineColor(kBlue);
  h_QE_Bkg->SetFillColor(kBlue);
  h_QE_Bkg->SetFillStyle(3444);

  MnvH1D* h_RES_Bkg = (MnvH1D*)bkg_file->Get((TString)name+"_RES");
  h_RES_Bkg->SetLineColor(kRed);
  h_RES_Bkg->SetFillColor(kRed);
  h_RES_Bkg->SetFillStyle(3444);

  MnvH1D* h_DIS_Bkg = (MnvH1D*)bkg_file->Get((TString)name+"_DIS");
  h_DIS_Bkg->SetLineColor(kGreen);
  h_DIS_Bkg->SetFillColor(kGreen);
  h_DIS_Bkg->SetFillStyle(3444);

  MnvH1D* h_2p2h_Bkg = (MnvH1D*)bkg_file->Get((TString)name+"_2p2h");
  h_2p2h_Bkg->SetLineColor(kBlack);
  h_2p2h_Bkg->SetFillColor(kBlack);
  h_2p2h_Bkg->SetFillStyle(3444);

  MnvH1D* h_Other_Bkg = (MnvH1D*)bkg_file->Get((TString)name+"_Other");
  h_Other_Bkg->SetLineColor(kGray);
  h_Other_Bkg->SetFillColor(kGray);
  h_Other_Bkg->SetFillStyle(3444);


  THStack* h = new THStack();
  h->Add(h_Other_Bkg);
  h->Add(h_2p2h_Bkg);
  h->Add(h_DIS_Bkg);
  h->Add(h_RES_Bkg);
  h->Add(h_QE_Bkg);

  h->Add(h_Other_Sig);
  h->Add(h_2p2h_Sig);
  h->Add(h_DIS_Sig);
  h->Add(h_RES_Sig);
  h->Add(h_QE_Sig);

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

  leg->AddEntry(h_QE_Sig,"Sig. + QE");
  leg->AddEntry(h_RES_Sig,"Sig. + RES");
  leg->AddEntry(h_DIS_Sig,"Sig. + DIS");
  leg->AddEntry(h_2p2h_Sig,"Sig. + 2p2h");
  leg->AddEntry(h_Other_Sig,"Sig. + Other");

  leg->AddEntry(h_QE_Bkg,"Bkg. + QE");
  leg->AddEntry(h_RES_Bkg,"Bkg. + RES");
  leg->AddEntry(h_DIS_Bkg,"Bkg. + DIS");
  leg->AddEntry(h_2p2h_Bkg,"Bkg. + 2p2h");
  leg->AddEntry(h_Other_Bkg,"Bkg. + Other");

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
  if (argc != 5) {
    cout << "Check usage..." << endl;
    return 2;
  }

  string sigName=string(argv[1]);
  string bkgName=string(argv[2]);
  string outDir=string(argv[3]);
  TString sample=argv[4];

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
  string sigNameStub = sigName;
  string bkgNameStub = bkgName;
  size_t pos=0;

  //cout << sigNameStub << endl;
  while ((pos = sigNameStub.find(slash)) != string::npos){
    //cout << sigNameStub << endl;
    token = sigNameStub.substr(0, pos);
    //cout << token << endl;
    sigNameStub.erase(0, pos+slash.length());
  }
  //cout << sigNameStub << endl;
  if ((pos=sigNameStub.find(rootExt)) == string::npos){
    cout << "Input need be .root file." << endl;
    return 4;
  }

  //cout << bkgNameStub << endl;
  while ((pos = bkgNameStub.find(slash)) != string::npos){
    //cout << bkgNameStub << endl;
    token = bkgNameStub.substr(0, pos);
    //cout << token << endl;
    bkgNameStub.erase(0, pos+slash.length());
  }
  //cout << bkgNameStub << endl;
  if ((pos=bkgNameStub.find(rootExt)) == string::npos){
    cout << "Input need be .root file." << endl;
    return 5;
  }

  cout << "Input Signal file name parsed to: " << sigNameStub << endl;

  TFile* sigFile = new TFile(sigName.c_str(),"READ");

  cout << "Input BKG file name parsed to: " << bkgNameStub << endl;

  TFile* bkgFile = new TFile(bkgName.c_str(),"READ");

  TList* keyList = sigFile->GetListOfKeys();
  if (!keyList){
    cout << "List of keys failed to get." << endl;
    return 6;
  }

  TIter next(keyList);
  TKey* key;
  while ( key = (TKey*)next() ){
    //cout << key->GetName() << endl;
    pos=0;
    string name=(string)key->GetName();
    if((pos = name.find("_QE")) == string::npos) continue;
    TCanvas* c1 = DrawToCanvas(name,sigFile,bkgFile,sample);
    name.erase(name.length()-3,name.length());
    c1->Print((TString)outDir+(TString)name+"_stacked_test.pdf");
    c1->Print((TString)outDir+(TString)name+"_stacked_test.png");
    delete c1;
  }

  //parsingTest("name",inFile);

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
