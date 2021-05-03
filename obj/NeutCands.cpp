//File: NeutCands.h
//Info: Neutron Candidate(s) Classes/NameSpace. Currently under development.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "NeutCands.h"

namespace NeutronCandidates{

  intBranchMap GetBranchIntMap(){
    return {{"SetID",{"MasterAnaDev_BlobID"}},{"3D",{"MasterAnaDev_BlobIs3D"}},{"SetMCPID",{"MasterAnaDev_BlobMCPID"}},{"SetTopMCPID",{"MasterAnaDev_BlobTopMCPID"}},{"SetMCParentTrackID",{"MasterAnaDev_BlobMCParentTrackID"}},};
  }

  doubleBranchMap GetBranchDoubleMap(){
    return {{"SetTotE",{"MasterAnaDev_BlobTotalE"}},{"SetBegPos",{"MasterAnaDev_BlobBegX","MasterAnaDev_BlobBegY","MasterAnaDev_BlobBegZ"}},{"SetEndPos",{"MasterAnaDev_BlobEndX","MasterAnaDev_BlobEndY","MasterAnaDev_BlobEndZ"}}, };
  }
  
  NeutCand::NeutCand(){
    this->init();
  }

  NeutCand::NeutCand(intCandData candIntData, doubleCandData candDoubleData, TVector3 vtx){
    this->init();
    this->SetEvtVtx(vtx);
    for(const auto& function: candIntData){
      if (function.first=="SetID") this->SetID(function.second);
      else if (function.first=="SetIs3D") this->SetIs3D(function.second);
      else if (function.first=="SetMCPID") this->SetMCPID(function.second);
      else if (function.first=="SetTopMCPID") this->SetTopMCPID(function.second);
      else continue;
    }
    for(const auto& function:candDoubleData){
      if (function.first=="SetTotE") this->SetTotalE(function.second);
      else if (function.first=="SetBegPos") this->SetBegPos(function.second);
      else if (function.first=="SetEndPos") this->SetEndPos(function.second);
      else continue;
    }
  }
  
  void NeutCand::init(){
    TVector3 tmp;
    tmp.SetXYZ(0.0,0.0,0.0);
    fID = -1;
    fIs3D = -999;
    fMCPID = -999;
    fTopMCPID = -999;
    fMCParentTrackID = -999;
    fTotE = -999.0;
    fAngleToFP = -999.0;
    fEvtVtx=tmp;
    fBegPos=tmp;
    fDirection=tmp;
    fFlightPath=tmp;
    tmp.~TVector3();
  }
  
  NeutCands::NeutCands(){
    this->init();
  }
  
  NeutCands::NeutCands(std::vector<NeutCand> cands){
    this->init();
    std::map<int, NeutCand> candsMap;
    double maxE = -1.0;
    for(int i_cand=0;i_cand < (int)cands.size();++i_cand){
      candsMap[cands.at(i_cand).GetID()]=cands.at(i_cand);
      if (cands.at(i_cand).GetTotalE() > maxE){
	fIDmaxE = cands.at(i_cand).GetID();
	fCandMaxE = cands.at(i_cand);
	maxE=cands.at(i_cand).GetTotalE();
      }
      this->SetCands(candsMap);
    }
  }
  
  NeutCands::NeutCands(std::map<int,intCandData> candsDataInt, std::map<int,doubleCandData> candDataDouble){
    this->init();
  }
  
  void NeutCands::init(){
    fCands = {};
    fNCands = fCands.size();
    fIDmaxE = -1;
    fCandMaxE = NeutCand();
  }

}
