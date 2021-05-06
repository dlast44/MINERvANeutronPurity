//File: EventLoop.cxx
//Info: This is a script to run a loop over all events in a single nTuple file and perform some plotting. Will eventually exist as the basis for the loops over events in analysis.
//
//Usage: EventLoop.cxx <MasterAnaDev_NTuple_list/single_file> <0=MC/1=PC> <output_directory> <tag_for_naming_files> optional: <n_event>
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

bool PassesFVCuts(CVUniverse& univ){
  vector<double> vtx = univ.GetVtx();
  double side = 850.0*2.0/sqrt(3.0);
  if (vtx[2] > targetBoundary || vtx[2] < 4500.0 ) return false;
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
  double leadingAngle=leadingBlob.GetAngleToFP();
  TVector3 leadingPos=leadingBlob.GetBegPos();
  TVector3 leadingFP=leadingBlob.GetFlightPath();
  TVector3 muonMom(univ.GetMuon4V().Z(),univ.GetMuon4V().Y(),univ.GetMuon4V().Z());
  double Q2GeV = univ.GetQ2QEPickledGeV();
  double MnGeV = 0.939566; //Should check what others might use, but this will be close enough for now...
  if (leadingFP.Mag()==0 || muonMom.Mag()==0) return 0;
  else{
    if((leading3D==1) && (leadingFP.Angle(muonMom) > 0.261799388)){
      if (leadingPos.Z() >= targetBoundary){
	if (leadingAngle > 0.2 && leadingAngle < 0.7) return 4;
	return 2;
      }
      if (leadingAngle > 0.2 && leadingAngle < 0.7) return 3;
      return 1;
    }
    return 0;
  }
}

bool PassesCuts(CVUniverse& univ, int isPC){
  //cout << "HELLO" << endl;
  return
    PassesFVCuts(univ) &&
    PassesCleanCCAntiNuCuts(univ, isPC) &&
    PassesTejinCCQECuts(univ);
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
  if (argc < 5 || argc > 6) {
    cout << "Check usage..." << endl;
    return 1;
  }

  string playlist=string(argv[1]);
  int isPC=atoi(argv[2]);
  string outDir=string(argv[3]);
  string tag=string(argv[4]);
  int nEntries=0;

  if (argc == 6){
    nEntries=atoi(argv[5]);
  }

  if (PathExists(outDir)){
    cout << "Thank you for choosing a path for output files that exists." << endl;
  }
  else{
    cout << "Output directory doesn't exist. Exiting" << endl;
    return 2;
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

  PlotUtils::HistWrapper<CVUniverse> hw_tracker_primary_parent_Tejin_TrackerONLY("hw_tracker_primary_parent_Tejin_TrackerONLY","Primary Particle Matched To Blob (CCQE, Recoil, Blob, Tracker Blob)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_tracker_primary_parent_Tejin("hw_tracker_primary_parent_Tejin","Primary Particle Matched To Blob (CCQE, Recoil, Blob)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_tracker_primary_parent_David_TrackerONLY("hw_tracker_primary_parent_David_TrackerONLY","Primary Particle Matched To Blob (CCQE, Recoil, Blob, Tracker Blob)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_tracker_primary_parent_David("hw_tracker_primary_parent_David","Primary Particle Matched To Blob (CCQE, Recoil, Blob)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_tracker_primary_parent_Recoil("hw_tracker_primary_parent_Recoil","Primary Particle Matched To Blob (CCQE, Recoil)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_tracker_primary_parent_CCQE("hw_tracker_primary_parent_CCQE","Primary Particle Matched To Blob (CCQE)",10,0,10,error_bands);

  PlotUtils::HistWrapper<CVUniverse> hw_target_primary_parent_Tejin_TrackerONLY("hw_target_primary_parent_Tejin_TrackerONLY","Primary Particle Matched To Blob (CCQE, Recoil, Blob, Tracker Blob)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_target_primary_parent_Tejin("hw_target_primary_parent_Tejin","Primary Particle Matched To Blob (CCQE, Recoil, Blob)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_target_primary_parent_David_TrackerONLY("hw_target_primary_parent_David_TrackerONLY","Primary Particle Matched To Blob (CCQE, Recoil, Blob, Tracker Blob)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_target_primary_parent_David("hw_target_primary_parent_David","Primary Particle Matched To Blob (CCQE, Recoil, Blob)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_target_primary_parent_Recoil("hw_target_primary_parent_Recoil","Primary Particle Matched To Blob (CCQE, Recoil)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_target_primary_parent_CCQE("hw_target_primary_parent_CCQE","Primary Particle Matched To Blob (CCQE)",10,0,10,error_bands);

  PlotUtils::HistWrapper<CVUniverse> hw_ALL_primary_parent_Tejin("hw_ALL_primary_parent_Tejin","Primary Particle Matched To Blob (CCQE, Recoil, Blob)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_ALL_primary_parent_David("hw_ALL_primary_parent_David","Primary Particle Matched To Blob (CCQE, Recoil, Blob)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_ALL_primary_parent_Tejin_TrackerONLY("hw_ALL_primary_parent_Tejin_TrackerONLY","Primary Particle Matched To Blob (CCQE, Recoil, Blob, Tracker Blob)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_ALL_primary_parent_David_TrackerONLY("hw_ALL_primary_parent_David_TrackerONLY","Primary Particle Matched To Blob (CCQE, Recoil, Blob)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_ALL_primary_parent_CCQE("hw_ALL_primary_parent_CCQE","Primary Particle Matched To Blob (CCQE)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_ALL_primary_parent_Recoil("hw_ALL_primary_parent_Recoil","Primary Particle Matched To Blob (CCQE, Recoil)",10,0,10,error_bands);

  if(!nEntries) nEntries = chain->GetEntries();
  cout << "Processing " << nEntries << " events." << endl;
  for (int i=0; i<nEntries;++i){
    if (i%(nEntries/100)==0) cout << (100*i)/nEntries << "% finished." << endl;
    //if (i%(10000)==0) cout << i << " entries finished." << endl;
    for (auto band : error_bands){
      vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes){
	universe->SetEntry(i);
	universe->UpdateNeutCands();
	
	//Passes CCQE Cuts that matche Tejin's selection
	if (PassesCuts(*universe, isPC)){
	  
	  //Passes Tejin Recoil and Blob
	  if (PassesTejinRecoilCut(*universe, isPC)){
	    
	    int TejinBlobValue = PassesTejinBlobCuts(*universe);
	    //Passes Tejin Recoil and Blob
	    if (TejinBlobValue){
	      NeutronCandidates::NeutCands cands = universe->GetCurrentNeutCands();
	      for (auto& cand: cands.GetCandidates()){
		int PID = cand.second.GetMCPID();
		int TopPID = cand.second.GetTopMCPID();
		int PTrackID = cand.second.GetMCParentTrackID();
		double candZ = cand.second.GetBegPos().Z();
		if (PTrackID==0 && !isPC){
		  //Additional Requirement of the Chosen Blob being in the tracker only.
		  if (TejinBlobValue==2){
		    hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->Fill(PDGbins[PID]);
		    if (candZ > targetBoundary){
		      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->Fill(PDGbins[PID]);
		    }
		    else {
		      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->Fill(PDGbins[PID]);
		    }
		  }
		  else if(TejinBlobValue==4){
		    hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->Fill(PDGbins[PID]);
		    hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->Fill(PDGbins[PID]);
		    hw_ALL_primary_parent_David.univHist(universe)->Fill(PDGbins[PID]);
		    if (candZ > targetBoundary){
		      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->Fill(PDGbins[PID]);
		      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->Fill(PDGbins[PID]);
		      hw_tracker_primary_parent_David.univHist(universe)->Fill(PDGbins[PID]);
		    }
		    else {
		      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->Fill(PDGbins[PID]);
		      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->Fill(PDGbins[PID]);
		      hw_target_primary_parent_David.univHist(universe)->Fill(PDGbins[PID]);
		    }
		  }
		  if (TejinBlobValue == 3) hw_ALL_primary_parent_David.univHist(universe)->Fill(PDGbins[PID]);
		  hw_ALL_primary_parent_Tejin.univHist(universe)->Fill(PDGbins[PID]);
		  hw_ALL_primary_parent_Recoil.univHist(universe)->Fill(PDGbins[PID]);
		  hw_ALL_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[PID]);
		  if (candZ > targetBoundary){
		    if (TejinBlobValue == 3) hw_tracker_primary_parent_David.univHist(universe)->Fill(PDGbins[PID]);
		    hw_tracker_primary_parent_Tejin.univHist(universe)->Fill(PDGbins[PID]);
		    hw_tracker_primary_parent_Recoil.univHist(universe)->Fill(PDGbins[PID]);
		    hw_tracker_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[PID]);
		  }
		  else{
		    if (TejinBlobValue == 3) hw_target_primary_parent_David.univHist(universe)->Fill(PDGbins[PID]); 
		    hw_target_primary_parent_Tejin.univHist(universe)->Fill(PDGbins[PID]);
		    hw_target_primary_parent_Recoil.univHist(universe)->Fill(PDGbins[PID]);
		    hw_target_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[PID]);
		  }
		}
		else{
		  //Additional Requirement of the Chosen Blob being in the tracker only.
		  if (TejinBlobValue==2){
		    hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->Fill(PDGbins[TopPID]);
		    if (candZ > targetBoundary){
		      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->Fill(PDGbins[TopPID]);
		    }		  
		    else {
		      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->Fill(PDGbins[TopPID]);
		    }
		  }
		  else if (TejinBlobValue==4){
		    hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->Fill(PDGbins[TopPID]);
		    hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->Fill(PDGbins[TopPID]);
		    hw_ALL_primary_parent_David.univHist(universe)->Fill(PDGbins[TopPID]);
		    if (candZ > targetBoundary){
		      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->Fill(PDGbins[TopPID]);
		      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->Fill(PDGbins[TopPID]);
		      hw_tracker_primary_parent_David.univHist(universe)->Fill(PDGbins[TopPID]);
		    }		  
		    else {
		      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->Fill(PDGbins[TopPID]);
		      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->Fill(PDGbins[TopPID]);
		      hw_target_primary_parent_David.univHist(universe)->Fill(PDGbins[TopPID]);
		    }
		  }
		  if (TejinBlobValue==3) hw_ALL_primary_parent_David.univHist(universe)->Fill(PDGbins[TopPID]);
		  hw_ALL_primary_parent_Tejin.univHist(universe)->Fill(PDGbins[TopPID]);
		  hw_ALL_primary_parent_Recoil.univHist(universe)->Fill(PDGbins[TopPID]);
		  hw_ALL_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[TopPID]);
		  if (candZ > targetBoundary){
		    if (TejinBlobValue==3) hw_tracker_primary_parent_David.univHist(universe)->Fill(PDGbins[TopPID]);
		    hw_tracker_primary_parent_Tejin.univHist(universe)->Fill(PDGbins[TopPID]);
		    hw_tracker_primary_parent_Recoil.univHist(universe)->Fill(PDGbins[TopPID]);
		    hw_tracker_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[TopPID]);
		  }
		  else{
		    if(TejinBlobValue==3) hw_target_primary_parent_David.univHist(universe)->Fill(PDGbins[TopPID]);
		    hw_target_primary_parent_Tejin.univHist(universe)->Fill(PDGbins[TopPID]);
		    hw_target_primary_parent_Recoil.univHist(universe)->Fill(PDGbins[TopPID]);
		    hw_target_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[TopPID]);
		  }
		}
	      }
	    }
	  
	    //Passes Tejin Recoil Not Blob
	    else {
	      NeutronCandidates::NeutCands cands = universe->GetCurrentNeutCands();
	      for (auto& cand: cands.GetCandidates()){
		//cout << "GOOD" << endl;	      
		int PID = cand.second.GetMCPID();
		int TopPID = cand.second.GetTopMCPID();
		int PTrackID = cand.second.GetMCParentTrackID();
		double candZ = cand.second.GetBegPos().Z();
		//if (cand.second.GetIs3D()==1) cout << "BlobIs3D" << endl;
		if (PTrackID==0 && !isPC){
		  hw_ALL_primary_parent_Recoil.univHist(universe)->Fill(PDGbins[PID]);
		  hw_ALL_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[PID]);
		  if (candZ > targetBoundary){
		    hw_tracker_primary_parent_Recoil.univHist(universe)->Fill(PDGbins[PID]);
		    hw_tracker_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[PID]);
		  }
		  else{
		    hw_target_primary_parent_Recoil.univHist(universe)->Fill(PDGbins[PID]);
		    hw_target_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[PID]);
		  }
		}
		else{
		  hw_ALL_primary_parent_Recoil.univHist(universe)->Fill(PDGbins[TopPID]);
		  hw_ALL_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[TopPID]);
		  if (candZ > targetBoundary){
		    hw_tracker_primary_parent_Recoil.univHist(universe)->Fill(PDGbins[TopPID]);
		    hw_tracker_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[TopPID]);
		  }
		  else{
		    hw_target_primary_parent_Recoil.univHist(universe)->Fill(PDGbins[TopPID]);
		    hw_target_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[TopPID]);
		  }
		}
	      }
	    }
	  }

	  //Fails Tejin Recoil. I'm not going to treat the Tejin Blob Cut as special/independent of this recoil cut.
	  else {
	    NeutronCandidates::NeutCands cands = universe->GetCurrentNeutCands();
	    for (auto& cand: cands.GetCandidates()){
	      //cout << "GOOD" << endl;	      
	      int PID = cand.second.GetMCPID();
	      int TopPID = cand.second.GetTopMCPID();
	      int PTrackID = cand.second.GetMCParentTrackID();
	      double candZ = cand.second.GetBegPos().Z();
	      //if (cand.second.GetIs3D()==1) cout << "BlobIs3D" << endl;
	      if (PTrackID==0 && !isPC){
		hw_ALL_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[PID]);
		if (candZ > targetBoundary) hw_tracker_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[PID]);
		else hw_target_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[PID]);
	      }
	      else{
		hw_ALL_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[TopPID]);
		if (candZ > targetBoundary) hw_tracker_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[TopPID]);
		else hw_target_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[TopPID]);
	      }
	    }
	  }
	}
      }
    }
  }

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->cd();
  for (auto band : error_bands){
    int i=0;
    vector<CVUniverse*> error_band_universes = band.second;
    for (auto universe : error_band_universes){
      ++i;

      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_Tejin.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_tracker_primary_parent_Tejin_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_tracker_primary_parent_Tejin_TrackerONLY_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_David.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_David.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_tracker_primary_parent_David.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_David.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_David.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_tracker_primary_parent_David_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_tracker_primary_parent_David_TrackerONLY_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_Recoil.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_tracker_primary_parent_Recoil_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");


      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_CCQE.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_tracker_primary_parent_CCQE_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+TString(to_string(nEntries))+"_Events.pdf");

      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_target_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_Tejin.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_target_primary_parent_Tejin_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_target_primary_parent_Tejin_TrackerONLY_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_David.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_David.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_target_primary_parent_David.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_David.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_David.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_target_primary_parent_David_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_target_primary_parent_David_TrackerONLY_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_target_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_Recoil.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_target_primary_parent_Recoil_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_target_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_CCQE.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_target_primary_parent_CCQE_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+TString(to_string(nEntries))+"_Events.pdf");

      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_Tejin.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_ALL_primary_parent_Tejin_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_Tejin_TrackerONLY.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_ALL_primary_parent_Tejin_TrackerONLY_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_David.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_David.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_ALL_primary_parent_David.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_David.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_David.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_ALL_primary_parent_David_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_David_TrackerONLY.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_ALL_primary_parent_David_TrackerONLY_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_Recoil.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_ALL_primary_parent_Recoil_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+(TString)(to_string(nEntries))+"_Events.pdf");

      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(1,"Other");
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(2,"");
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(3,"n");           
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(4,"p");           
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(9,"e");
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetBinLabel(10,"#mu");
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetLabelSize(0.06);
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetXaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetTitle("Blobs");
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_CCQE.univHist(universe)->Draw();
      c1->Print((TString)(outDir)+"h_ALL_primary_parent_CCQE_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+TString(to_string(nEntries))+"_Events.pdf");

      //Drawing target vs. tracker on same plot instead of stacked.

      hw_target_primary_parent_CCQE.univHist(universe)->SetLineColor(kBlue);
      //hw_target_primary_parent_CCQE.univHist(universe)->Sumw2(kFALSE);
      hw_target_primary_parent_CCQE.univHist(universe)->Scale(1.0/(((double)hw_target_primary_parent_CCQE.univHist(universe)->Integral())));
      hw_target_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_target_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetRangeUser(0.0,0.5);
      hw_tracker_primary_parent_CCQE.univHist(universe)->SetLineColor(kRed);
      //hw_tracker_primary_parent_CCQE.univHist(universe)->Sumw2(kFALSE);
      hw_tracker_primary_parent_CCQE.univHist(universe)->Scale(1.0/(((double)hw_tracker_primary_parent_CCQE.univHist(universe)->Integral())));
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_tracker_primary_parent_CCQE.univHist(universe)->GetYaxis()->SetRangeUser(0.0,0.5);

      hw_target_primary_parent_Recoil.univHist(universe)->SetLineColor(kBlue);
      //hw_target_primary_parent_Recoil.univHist(universe)->Sumw2(kFALSE);
      hw_target_primary_parent_Recoil.univHist(universe)->Scale(1.0/(((double)hw_target_primary_parent_Recoil.univHist(universe)->Integral())));
      hw_target_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_target_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetRangeUser(0.0,0.5);
      hw_tracker_primary_parent_Recoil.univHist(universe)->SetLineColor(kRed);
      //hw_tracker_primary_parent_Recoil.univHist(universe)->Sumw2(kFALSE);
      hw_tracker_primary_parent_Recoil.univHist(universe)->Scale(1.0/(((double)hw_tracker_primary_parent_Recoil.univHist(universe)->Integral())));
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_tracker_primary_parent_Recoil.univHist(universe)->GetYaxis()->SetRangeUser(0.0,0.5);

      hw_target_primary_parent_Tejin.univHist(universe)->SetLineColor(kBlue);
      //hw_target_primary_parent_Tejin.univHist(universe)->Sumw2(kFALSE);
      hw_target_primary_parent_Tejin.univHist(universe)->Scale(1.0/(((double)hw_target_primary_parent_Tejin.univHist(universe)->Integral())));
      hw_target_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_target_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetRangeUser(0.0,0.5);
      hw_tracker_primary_parent_Tejin.univHist(universe)->SetLineColor(kRed);
      //hw_tracker_primary_parent_Tejin.univHist(universe)->Sumw2(kFALSE);
      hw_tracker_primary_parent_Tejin.univHist(universe)->Scale(1.0/(((double)hw_tracker_primary_parent_Tejin.univHist(universe)->Integral())));
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_tracker_primary_parent_Tejin.univHist(universe)->GetYaxis()->SetRangeUser(0.0,0.5);

      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->SetLineColor(kBlue);
      //hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->Sumw2(kFALSE);
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->Scale(1.0/(((double)hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->Integral())));
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetRangeUser(0.0,0.5);
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->SetLineColor(kRed);
      //hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->Sumw2(kFALSE);
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->Scale(1.0/(((double)hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->Integral())));
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->GetYaxis()->SetRangeUser(0.0,0.5);

      hw_target_primary_parent_David.univHist(universe)->SetLineColor(kBlue);
      //hw_target_primary_parent_David.univHist(universe)->Sumw2(kFALSE);
      hw_target_primary_parent_David.univHist(universe)->Scale(1.0/(((double)hw_target_primary_parent_David.univHist(universe)->Integral())));
      hw_target_primary_parent_David.univHist(universe)->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_target_primary_parent_David.univHist(universe)->GetYaxis()->SetRangeUser(0.0,0.5);
      hw_tracker_primary_parent_David.univHist(universe)->SetLineColor(kRed);
      //hw_tracker_primary_parent_David.univHist(universe)->Sumw2(kFALSE);
      hw_tracker_primary_parent_David.univHist(universe)->Scale(1.0/(((double)hw_tracker_primary_parent_David.univHist(universe)->Integral())));
      hw_tracker_primary_parent_David.univHist(universe)->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_tracker_primary_parent_David.univHist(universe)->GetYaxis()->SetRangeUser(0.0,0.5);

      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->SetLineColor(kBlue);
      //hw_target_primary_parent_David_TrackerONLY.univHist(universe)->Sumw2(kFALSE);
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->Scale(1.0/(((double)hw_target_primary_parent_David_TrackerONLY.univHist(universe)->Integral())));
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetRangeUser(0.0,0.5);
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->SetLineColor(kRed);
      //hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->Sumw2(kFALSE);
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->Scale(1.0/(((double)hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->Integral())));
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->GetYaxis()->SetRangeUser(0.0,0.5);

      TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
      legend->SetHeader("Blob Location");
      legend->AddEntry(hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe),"Tracker Region","L");
      legend->AddEntry(hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe),"Target Region","L");

      hw_tracker_primary_parent_CCQE.univHist(universe)->SetStats(kFALSE);
      hw_target_primary_parent_CCQE.univHist(universe)->SetStats(kFALSE);
      hw_tracker_primary_parent_CCQE.univHist(universe)->Draw();
      hw_target_primary_parent_CCQE.univHist(universe)->Draw("same");
      legend->Draw();
      c1->Print((TString)(outDir)+"h_both_regions_primary_parent_CCQE_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+TString(to_string(nEntries))+"_Events.pdf");

      hw_tracker_primary_parent_Recoil.univHist(universe)->SetStats(kFALSE);
      hw_target_primary_parent_Recoil.univHist(universe)->SetStats(kFALSE);
      hw_tracker_primary_parent_Recoil.univHist(universe)->Draw();
      hw_target_primary_parent_Recoil.univHist(universe)->Draw("same");
      legend->Draw();
      c1->Print((TString)(outDir)+"h_both_regions_primary_parent_Recoil_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+TString(to_string(nEntries))+"_Events.pdf");

      hw_tracker_primary_parent_Tejin.univHist(universe)->SetStats(kFALSE);
      hw_target_primary_parent_Tejin.univHist(universe)->SetStats(kFALSE);
      hw_tracker_primary_parent_Tejin.univHist(universe)->Draw();
      hw_target_primary_parent_Tejin.univHist(universe)->Draw("same");
      legend->Draw();
      c1->Print((TString)(outDir)+"h_both_regions_primary_parent_Tejin_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+TString(to_string(nEntries))+"_Events.pdf");

      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->SetStats(kFALSE);
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->SetStats(kFALSE);
      hw_tracker_primary_parent_Tejin_TrackerONLY.univHist(universe)->Draw();
      hw_target_primary_parent_Tejin_TrackerONLY.univHist(universe)->Draw("same");
      legend->Draw();
      c1->Print((TString)(outDir)+"h_both_regions_primary_parent_Tejin_TrackerONLY_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+TString(to_string(nEntries))+"_Events.pdf");

      hw_tracker_primary_parent_David.univHist(universe)->SetStats(kFALSE);
      hw_target_primary_parent_David.univHist(universe)->SetStats(kFALSE);
      hw_tracker_primary_parent_David.univHist(universe)->Draw();
      hw_target_primary_parent_David.univHist(universe)->Draw("same");
      legend->Draw();
      c1->Print((TString)(outDir)+"h_both_regions_primary_parent_David_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+TString(to_string(nEntries))+"_Events.pdf");

      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->SetStats(kFALSE);
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->SetStats(kFALSE);
      hw_tracker_primary_parent_David_TrackerONLY.univHist(universe)->Draw();
      hw_target_primary_parent_David_TrackerONLY.univHist(universe)->Draw("same");
      legend->Draw();
      c1->Print((TString)(outDir)+"h_both_regions_primary_parent_David_TrackerONLY_"+(TString)(universe->ShortName())+(TString)(to_string(i))+"_"+TString(playlistStub)+"_"+TString(tag)+"_"+TString(to_string(nEntries))+"_Events.pdf");

    }
  }

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;

}
