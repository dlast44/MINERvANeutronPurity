//File: CVUniverse.h
//Info: Central Value Universe Class. Currently under development.
//      Utilized in following New Systematics Framework (MAT) approach.
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

#include "PlotUtils/MinervaUniverse.h"
#include "obj/NeutCands.h"
#include "TVector3.h"

class CVUniverse: public PlotUtils::MinervaUniverse {
 public:
  //CTOR
 CVUniverse(typename PlotUtils::MinervaUniverse::config_t chw, const double nsigma=0): PlotUtils::MinervaUniverse(chw, nsigma) {}

  //DTOR
  virtual ~CVUniverse() = default;

  //Shared (pre-defined, common, etc.) systematics components
  #include "PlotUtils/SystCalcs/WeightFunctions.h"
  #include "PlotUtils/SystCalcs/MuonFunctions.h"
  #include "PlotUtils/SystCalcs/TruthFunctions.h"

  //Useful naming grab based on the inherent object in the class iself that should work *crosses fingers*
  //virtual std::string GetAnaToolName() const { return (std::string)m_chw->GetName(); }

  //Initial Reco Branches to investigate
  virtual int GetNTracks() const { return GetInt("multiplicity"); };
  virtual int GetNNeutBlobs() const { return GetInt("MasterAnaDev_BlobIs3D_sz"); };
  virtual std::vector<double> GetEMBlobStartZVec() const { return GetVec<double>("nonvtx_iso_blobs_start_position_z_in_prong"); };
  virtual std::vector<int> GetEMBlobNHitsVec() const { return GetVec<int>("nonvtx_iso_blobs_n_hits_in_prong"); };
  virtual std::vector<double> GetEMBlobEnergyVec() const { return GetVec<double>("nonvtx_iso_blobs_energy_in_prong"); };
  virtual std::vector<double> GetEMNBlobsTotalEnergyTotalNHits(double shift = 0) const {
    std::vector<double> info;
    double nBlobs = 0;
    double totalE = shift;
    double nHits = 0;
    std::vector<double> StartZVec = GetEMBlobStartZVec();
    std::vector<double> EnergyVec = GetEMBlobEnergyVec();
    std::vector<int> NHitsVec = GetEMBlobNHitsVec();
    for (unsigned int i=0; i<StartZVec.size(); ++i){
      if (StartZVec.at(i) > 4750.0){
	nBlobs+=1.0;
	totalE+=EnergyVec.at(i);
	nHits+=(double)NHitsVec.at(i);
      }
    }
    info.push_back(nBlobs);
    info.push_back(totalE);
    info.push_back(nHits);
    return info;
  };

  virtual std::vector<double> GetVtx() const { return GetVec<double>("vtx"); };

  virtual int GetNImprovedMichel() const { return GetInt("improved_michel_vertex_type_sz"); };

  virtual int GetNDeadDiscriminatorsUpstreamMuon() const { return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj"); };
  virtual int GetIsMinosMatchTrack() const { return GetInt("isMinosMatchTrack"); };
  virtual int GetIsMinosMatchTrackOLD() const { return GetInt("muon_is_minos_match_track"); };
  virtual int GetIsMinosMatchStub() const { return GetInt("isMinosMatchStub"); };
  virtual int GetIsMinosMatchStubOLD() const { return GetInt("muon_is_minos_match_stub"); };
  virtual int GetNuHelicity() const { return GetInt("MasterAnaDev_nuHelicity"); };

  //Neutron Candidate Business
  virtual int GetNNeutBlobs() const { return GetInt("MasterAnaDev_BlobID_sz"); };

  virtual NeutronCandidates::NeutCand GetNeutCand(int index) const{
    std::vector<double> vtx = GetVtx();
    TVector3 EvtVtx;
    EvtVtx.setXYZ(vtx.at(0),vtx.at(1),vtx.at(2));
    NeutronCandidates::intCandData intData;
    NeutronCandidates::doubleCandData doubleData;
    for (const auto& intMember: NeutronCandidates::IntBranchList){
      intData[intMember.first]={}
      for (const auto& branchName: intMember.second){
	intData[intMember.first].push_back(GetVecElemInt("MasterAnaDev"+branchName),index);
      }
    }
    for (const auto& doubleMember: NeutronCandidates::DoubleBranchList){
      doubleData[doubleMember.first]={}
      for (const auto& branchName: doubleMember.second){
	doubleData[doubleMember.first].push_back(GetVecElem("MasterAnaDev"+branchName),index);
      }
    }
    return NeutronCandidates::NeutCand(intData,doubleData,EvtVtx);
  };
  
  virtual NeutronCandidates::NeutCands GetNeutCands() const{
    std::vector<NeutronCandidates::NeutCand> cands;
    int nBlobs = GetNNeutBlobs();
    for(int neutBlobIndex=0; neutBlobIndex < nBlobs; ++neutBlobIndex){
      cands.push_back(GetNeutCand(neutBlobIndex));
    }
    NeutronCandidates::NeutCands EvtCands(cands);
    return EvtCands;
  };

  virtual void UpdateNeutCands() const{
    fNeutCands = GetNeutCands();
  };

  virtual bool CheckNeutBlobs() const{
    return fNeutCands.GetCandidates.size()==GetNNeutBlobs();
  }

 private:
  NeutronCandidates::NeutCands fNeutCands;
  
};

#endif
