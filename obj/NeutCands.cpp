//File: NeutCands.h
//Info: Neutron Candidate(s) Classes/NameSpace. Currently under development.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "NeutCands.h"
#include "TVector3.h"

using namespace NeutronCandidates;

NeutCand(){}

NeutCand(intCandData candIntData, doubleCandData candDoubleData, TVector3 vtx){
  this.SetEvtVtx(vtx);
  for(const auto& function:candIntData){
    if (function.first=="SetID") this.SetID(function.second);
    else if (function.first=="SetIs3D") this.SetIs3D(function.second);
    else continue;
  }
  for(const auto& function:candIntData){
    if (function.first=="SetTotE") this.SetTotE(function.second);
    else if (function.first=="SetBegPos") this.SetBegPos(function.second);
    else if (function.first=="SetEndPos") this.SetEndPos(function.second);
    else continue;
  }
}

NeutCands(){}

NeutCands(std::vector<NeutCand> cands){
  std::map<>
}

NeutCands(std::map<int,intCandData> candsDataInt, std::map<int,doubleCandData> candDataDouble){
  
}

}
#endif
