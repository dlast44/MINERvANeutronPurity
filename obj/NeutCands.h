//File: NeutCands.h
//Info: Neutron Candidate Class. Currently under development.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef NEUTCANDS_H
#define NEUTCANDS_H

#include "TVector3.h"
#include "stdlib.h"
#include <string>
#include <vector>
#include <map>

namespace NeutronCandidates{
  typedef std::map<std::string, std::vector<const char*>> intBranchMap;
  typedef std::map<std::string, std::vector<const char*>> doubleBranchMap;
  typedef std::map<std::string, std::vector<int>> intCandData;
  typedef std::map<std::string, std::vector<double>> doubleCandData;

  intBranchMap GetBranchIntMap();
  doubleBranchMap GetBranchDoubleMap();

  class NeutCand{
  private:
    //Currently only coding in the members that I actively use in MnvTgtNeutrons/particleCannon/nonMAT/interactiveMacros/Basic_Cuts_Try.cc
    
    int fID;
    int fIs3D;
    int fMCPID;
    int fTopMCPID;
    int fMCParentTrackID;
    double fTotE;
    double fAngleToFP;
    TVector3 fEvtVtx;
    TVector3 fBegPos;
    TVector3 fEndPos;
    TVector3 fDirection;
    TVector3 fFlightPath;

    void init();

  public:
    //CTOR
    NeutCand();
    NeutCand(NeutronCandidates::intCandData candIntData, NeutronCandidates::doubleCandData candDoubleData, TVector3 vtx);

    int GetID(){ return fID; };
    int GetIs3D(){ return fIs3D; };
    int GetMCPID(){ return fMCPID; };
    int GetTopMCPID(){ return fTopMCPID; };
    int GetMCParentTrackID(){ return fMCParentTrackID; };
    double GetTotalE(){ return fTotE; };
    double GetAngleToFP(){ return fAngleToFP; };
    TVector3 GetBegPos(){ return fBegPos; };
    TVector3 GetFlightPath(){ return fFlightPath; };
    TVector3 GetDirection(){ return fDirection; };
    TVector3 GetEvtVtx(){ return fEvtVtx; };

    void SetID(std::vector<int> ID){ fID=ID.at(0); };
    void SetIs3D(std::vector<int> is3D){ fIs3D=is3D.at(0); };
    void SetMCPID(std::vector<int> MCPID){ fMCPID=MCPID.at(0); };
    void SetTopMCPID(std::vector<int> TopPID){ fTopMCPID=TopPID.at(0); };
    void SetMCParentTrackID(std::vector<int> ParentID){ fTopMCPID=ParentID.at(0); };
    void SetTotalE(std::vector<double> TotE){ fTotE=TotE.at(0); };
    void SetEvtVtx(TVector3 EvtVtx){ fEvtVtx=EvtVtx; };
    //Move MULTI-LINE DEFINITIONS TO CPP...???
    void SetBegPos(std::vector<double> BegPos){
      fBegPos.SetXYZ(BegPos.at(0),BegPos.at(1),BegPos.at(2));
      fDirection = fBegPos-fEndPos;
      fFlightPath = fBegPos-fEvtVtx;
      if (fFlightPath.Mag() > 0 && fDirection.Mag() > 0){
	fAngleToFP = fFlightPath.Angle(fDirection);
      }
      else {
	fAngleToFP = -9999.0;
      }
    };
    void SetEndPos(std::vector<double> EndPos){
      fEndPos.SetXYZ(EndPos.at(0),EndPos.at(1),EndPos.at(2));
      fDirection = fBegPos-fEndPos;
      fAngleToFP = fFlightPath.Angle(fDirection);
      if (fFlightPath.Mag() > 0 && fDirection.Mag() > 0){
	fAngleToFP = fFlightPath.Angle(fDirection);
      }
      else {
	fAngleToFP = -9999.0;
      }
    };

    //DTOR
    virtual ~NeutCand() = default;
  };

  class NeutCands {
  private:
    int fNCands;
    int fIDmaxE;
    NeutCand fCandMaxE;
    std::map<int, NeutCand> fCands;

    void init();

  public:
    //CTORS
    NeutCands();
    NeutCands(std::vector<NeutCand> cands);
    NeutCands(std::map<int,intCandData> candsDataInt, std::map<int,doubleCandData> candDataDouble);

    //DTOR
    virtual ~NeutCands() = default;
    
    void SetCands(std::map<int, NeutCand> inCands){ 
      fCands=inCands; 
      fNCands=inCands.size();
    }

    int GetIDMaxE(){ return fIDmaxE; };
    int GetNCands(){ return fNCands; };
    NeutCand GetCandidate(int ID){ 
      if (fNCands == 0) return NeutCand();
      else return fCands[ID]; };
    NeutCand GetMaxCandidate(){ return fCandMaxE; };
    std::map<int, NeutCand> GetCandidates(){ return fCands; };
  };
}
#endif
