//File: NeutCands.h
//Info: Neutron Candidate Class. Currently under development.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef NEUTCANDS_H
#define NEUTCANDS_H

#include "TVector3.h"

namespace NeutronCandidates{
  typedef intBranchMap std::map<std::string, std::vector<std::string>>;
  typedef doubleBranchMap std::map<std::string, std::vector<std::string>>;
  typedef intCandData std::map<std::string, std::vector<int>>;
  typedef doubleCandData std::map<std::string, std::vector<double>>;

  intBranchMap IntBranchMap{{"SetID",{"_BlobID"}}, {"SetIs3D",{"_BlobIs3D"}}};
  doubleBranchMap DoubleBranchMap{{"SetTotE","_BlobTotalE"},{"SetBegPos",{"_BlobBegX","_BlobBegY","_BlobBegZ"}},{"SetEndPos",{"_BlobEndX","_BlobEndY","_BlobEndZ"}}};

  class NeutCand{
  public:
    //CTOR
    NeutCand();
    NeutCand(intCandData candIntData, doubleCandData candDoubleData, TVector3 vtx);

    void SetID(std::vector<int> ID){ fID=ID.at(0); };
    void SetIs3D(std::vector<int> is3D){ fIs3D=is3D.at(0); };
    void SetTotalE(std::vector<double> TotE){ fTotE=totE.at(0); };
    void SetEvtVtx(std::vector<double>EvtVtx){ fEvtVtx.SetXYZ(EvtVtx.at(0),EvtVtx.at(1),EvtVtx.at(2)); };
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
  private:
    //Currently only coding in the members that I actively use in MnvTgtNeutrons/particleCannon/nonMAT/interactiveMacros/Basic_Cuts_Try.cc
    
    int fID;
    int fIs3D;
    double fTotE;
    double fAngleToFP = -999.0;
    TVector3 fEvtVtx;
    TVector3 fBegPos;
    TVector3 fEndPos;
    TVector3 fDirection;
    TVector3 fFlightPath;
  };

  class NeutCands {
  public:
    //CTORS
    NeutCands();
    NeutCands(std::vector<NeutCand> cands);
    NeutCands(std::map<int,intCandData> candsDataInt, std::map<int,doubleCandData> candDataDouble);

    //DTOR
    virtual ~NeutCands() = default;

    void SetCands(std::map<int, NeutCand> inCands){ fCands=inCands; }

    int GetNCands(){ return fNCands; };
    NeutCand GetCandidate(int ID){ return fCands[i]};
    std::map<std::int, NeutCand> GetCandidates(){ return fCands; };
  private:
    int fNCands;
    int fIDmaxE = -1;
    std::map<int, NeutCand> fCands;
  };
}
#endif
