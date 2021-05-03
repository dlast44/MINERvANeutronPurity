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
  int leadingBlobID=universe.GetCurrentNeutCands().GetIDMaxE();
  //cout << "leading Blob ID " << leadingBlobID << endl;
  NeutronCandidates::NeutCand leadingBlob=universe.GetCurrentNeutCands().GetCandidate(leadingBlobID);

  //NeutronCandidates::NeutCand leadingBlob=universe.GetCurrentLeadingNeutCand();
  TVector3 leadingPos=leadingBlob.GetBegPos();
  TVector3 leadingFP=leadingBlob.GetFlightPath();
  TVector3 muonMom(universe.GetMuon4V().Z(),universe.GetMuon4V().Y(),universe.GetMuon4V().Z());
  if (leadingFP.Mag()==0 || muonMom.Mag()==0) return false;
  else{
    return
      (fabs(leadingPos.X()+999.0) > 0.5) &&
      (fabs(leadingPos.Y()+999.0) > 0.5) &&
      (fabs(leadingPos.Z()+999.0) > 0.5) &&
      (leadingFP.Angle(muonMom) > 0.261799388);
  }
}

bool PassesTejinBlobCutsTEST(CVUniverse& universe){
  NeutronCandidates::NeutCand leadingBlob=universe.GetCurrentLeadingNeutCand();
  TVector3 leadingPos=leadingBlob.GetBegPos();
  TVector3 leadingFP=leadingBlob.GetFlightPath();
  TVector3 muonMom(universe.GetMuon4V().Z(),universe.GetMuon4V().Y(),universe.GetMuon4V().Z());
  if (leadingFP.Mag()==0 || muonMom.Mag()==0) return false;
  else{
    return
      (fabs(leadingPos.X()+999.0) > 0.5) &&
      (fabs(leadingPos.Y()+999.0) > 0.5) &&
      (fabs(leadingPos.Z()+999.0) > 0.5) &&
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
  map< string, vector<CVUniverse*>> error_bands;
  error_bands[string("CV")].push_back(CV);

  //PlotUtils::HistWrapper<CVUniverse> hw_nBlobs("hw_nBlobs","Number of Neutron Blobs for this selection MAD 6D \"Full\" (04-28-2021 10:11 AM)",100,0.0,100.0,error_bands);
  //PlotUtils::HistWrapper<CVUniverse> hw_nBlobs_from_CandObj("hw_nBlobs_fromCandObj","Number of Neutron Blobs for this selection MAD 6D \"Full\" (04-28-2021 10:11 AM)",100,0.0,100.0,error_bands);

  PlotUtils::HistWrapper<CVUniverse> hw_primary_parent_Tejin("hw_primary_parent_Tejin","Primary Particle Matched To Blob (Tejin Selection w/ Blob Cut)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_primary_parent_CCQE("hw_primary_parent_CCQE","Primary Particle Matched To Blob (Tejin Selection w/o Blob Cut)",10,0,10,error_bands);
  PlotUtils::HistWrapper<CVUniverse> hw_primary_parent_Tejin_TEST("hw_primary_parent_Tejin_TEST","Primary Particle Matched To Blob (Tejin Selection w/ Blob Cut)",10,0,10,error_bands);

  int nEntries = chain->GetEntries();
  cout << "Processing " << nEntries << " events." << endl;
  for (int i=0; i<nEntries;++i){
    //if (i%(nEntries/100)==0) cout << (100*i)/nEntries << "% finished." << endl;
    for (auto band : error_bands){
      vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes){
	universe->SetEntry(i);
	universe->UpdateNeutCands();
	cout << "N Cands Before General Cuts A: " << universe->GetNNeutCands() << endl;
	cout << "N Cands Before General Cuts B: " << universe->GetCurrentNeutCands().GetNCands() << endl;
	cout << "N Cands Before General Cuts C: " << universe->GetCurrentNeutCands().GetCandidates().size() << endl;

	if (PassesCuts(*universe, isPC)){
	  //else{ This was here before... Checking that this has nothing to do with the changed behavior...
	  //There's weird buggy behavior...
	  NeutronCandidates::NeutCands cands = universe->GetCurrentNeutCands();
	  for (auto& cand: cands.GetCandidates()){
	    //cout << "GOOD" << endl;	      
	    int PID = cand.second.GetMCPID();
	    int TopPID = cand.second.GetTopMCPID();
	    int PTrackID = cand.second.GetMCParentTrackID();
	    //if (cand.second.GetIs3D()==1) cout << "BlobIs3D" << endl;
	    if (PTrackID==0 && !isPC){
	      hw_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[PID]);
	    }
	    else{
	      hw_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[TopPID]);
	    }
	  }

	cout << "N Cands Passed General Cuts, Before default blob, A: " << universe->GetNNeutCands() << endl;
	cout << "N Cands Passed General Cuts, Before default blob, B: " << universe->GetCurrentNeutCands().GetNCands() << endl;
	cout << "N Cands Passed General Cuts, Before default blob, C: " << universe->GetCurrentNeutCands().GetCandidates().size() << endl;

	  //cout << "Passes CCQE" << endl;
	  if (PassesTejinBlobCuts(*universe)){
	    //cout << "Passes Tejin" << endl;
	    NeutronCandidates::NeutCands cands = universe->GetCurrentNeutCands();
	    for (auto& cand: cands.GetCandidates()){
	      //cout << "GOOD" << endl;	      
	      int PID = cand.second.GetMCPID();
	      int TopPID = cand.second.GetTopMCPID();
	      int PTrackID = cand.second.GetMCParentTrackID();
	      if (PTrackID==0 && !isPC){
		hw_primary_parent_Tejin.univHist(universe)->Fill(PDGbins[PID]);
		//hw_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[PID]);
	      }
	      else{
		hw_primary_parent_Tejin.univHist(universe)->Fill(PDGbins[TopPID]);
		//hw_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[TopPID]);
	      }
	    }
	  }

	cout << "N Cands Passed General Cuts, After default blob, A: " << universe->GetNNeutCands() << endl;
	cout << "N Cands Passed General Cuts, After default blob, B: " << universe->GetCurrentNeutCands().GetNCands() << endl;
	cout << "N Cands Passed General Cuts, After default blob, C: " << universe->GetCurrentNeutCands().GetCandidates().size() << endl;

	  if (PassesTejinBlobCutsTEST(*universe)){
	    //cout << "Passes Tejin" << endl;
	    NeutronCandidates::NeutCands cands = universe->GetCurrentNeutCands();
	    for (auto& cand: cands.GetCandidates()){
	      //cout << "GOOD" << endl;	      
	      int PID = cand.second.GetMCPID();
	      int TopPID = cand.second.GetTopMCPID();
	      int PTrackID = cand.second.GetMCParentTrackID();
	      if (PTrackID==0 && !isPC){
		hw_primary_parent_Tejin_TEST.univHist(universe)->Fill(PDGbins[PID]);
		//hw_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[PID]);
	      }
	      else{
		hw_primary_parent_Tejin_TEST.univHist(universe)->Fill(PDGbins[TopPID]);
		//hw_primary_parent_CCQE.univHist(universe)->Fill(PDGbins[TopPID]);
	      }
	    }
	  }

	cout << "N Cands Passed General Cuts, After All blob, A: " << universe->GetNNeutCands() << endl;
	cout << "N Cands Passed General Cuts, After All blob, B: " << universe->GetCurrentNeutCands().GetNCands() << endl;
	cout << "N Cands Passed General Cuts, After All blob, C: " << universe->GetCurrentNeutCands().GetCandidates().size() << endl;

	  //hw_nBlobs.univHist(universe)->Fill(universe->GetNNeutBlobs());
	  //hw_nBlobs_from_CandObj.univHist(universe)->Fill(universe->GetNNeutCands());
	}
      }
    }
  }

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;

}
