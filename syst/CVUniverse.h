//File: CVUniverse.h
//Info: Central Value Universe Class. Currently under development.
//      Utilized in following New Systematics Framework (MAT) approach.
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

#include "PlotUtils/PhysicsVariables.h"
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
  #include "PlotUtils/SystCalcs/RecoilEnergyFunctions.h"

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

  virtual int GetHasInteractionVertex() const { return GetInt("has_interaction_vertex"); };

  virtual std::vector<double> GetVtx() const { return GetVec<double>("vtx"); };

  virtual int GetNFSPart() const { return GetInt("mc_nFSPart"); };

  virtual int GetNImprovedMichel() const { return GetInt("improved_michel_vertex_type_sz"); };

  virtual int GetNDeadDiscriminatorsUpstreamMuon() const { return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj"); };
  virtual int GetIsMinosMatchTrack() const { return GetInt("isMinosMatchTrack"); };
  virtual int GetIsMinosMatchTrackOLD() const { return GetInt("muon_is_minos_match_track"); };
  virtual int GetIsMinosMatchStub() const { return GetInt("isMinosMatchStub"); };
  virtual int GetIsMinosMatchStubOLD() const { return GetInt("muon_is_minos_match_stub"); };
  virtual int GetNuHelicity() const { return GetInt("MasterAnaDev_nuHelicity"); };

  double MeVGeV=0.001;

  virtual double GetCalRecoilEnergy() const{
    if (GetVec<double>("recoil_summed_energy").size()==0) return -999.0;
    return (GetVec<double>("recoil_summed_energy")[0]-GetDouble("recoil_energy_nonmuon_nonvtx100mm"));
  };

  virtual double GetNonCalRecoilEnergy() const{
    return 0.0;
  }
  
  virtual double GetEnuCCQEPickledGeV() const{ //RETURNS IN MeV^2
    int charge=-1; //hard-coded since I'm focused on anti-nu
    double enu=PlotUtils::nuEnergyCCQE( GetEmu(), GetPmu(), GetThetamu(), charge)*MeVGeV;
    return enu;
  };

  virtual double GetQ2QEPickledGeV() const{ //RETURNS IN MeV^2
    int charge=-1; //hard-coded since I'm focused on anti-nu
    if (GetEnuCCQEPickledGeV()<=0.0) return 0.0;
    else{
      double q2=PlotUtils::qSquaredCCQE( GetEmu(), GetPmu(), GetThetamu(), charge)*MeVGeV*MeVGeV;
      return q2;
    }
  }

  virtual double GetRecoilEnergyGeV() const{
    double recoilE = GetRecoilEnergy()*MeVGeV;
    return recoilE;
  }

  //Neutron Candidate Business

  virtual NeutronCandidates::NeutCand GetNeutCand(int index){
    std::vector<double> vtx = GetVtx();
    TVector3 EvtVtx;
    EvtVtx.SetXYZ(vtx.at(0),vtx.at(1),vtx.at(2));
    NeutronCandidates::intCandData intData;
    NeutronCandidates::doubleCandData doubleData;
    for (const auto& intMember: NeutronCandidates::GetBranchIntMap()){
      intData[intMember.first]={};
      for (const auto& branchName: intMember.second){
	intData[intMember.first].push_back(GetVecElemInt(branchName,index));
      }
    }
    for (const auto& doubleMember: NeutronCandidates::GetBranchDoubleMap()){
      doubleData[doubleMember.first]={};
      for (const auto& branchName: doubleMember.second){
	doubleData[doubleMember.first].push_back(GetVecElem(branchName,index));
      }
    }
    return NeutronCandidates::NeutCand(intData,doubleData,EvtVtx);
  };
  
  virtual NeutronCandidates::NeutCands GetNeutCands(){
    std::vector<NeutronCandidates::NeutCand> cands;
    int nBlobs = GetNNeutBlobs();
    for(int neutBlobIndex=0; neutBlobIndex < nBlobs; ++neutBlobIndex){
      cands.push_back(GetNeutCand(neutBlobIndex));
    }
    NeutronCandidates::NeutCands EvtCands(cands);
    return EvtCands;
  };

  virtual void UpdateNeutCands(){
    NeutronCandidates::NeutCands candsIn = GetNeutCands();
    fNeutCands = candsIn;
    fNNeutCands = candsIn.GetNCands();
  };

  NeutronCandidates::NeutCand GetCurrentNeutCand(int index) { return fNeutCands.GetCandidate(index); };

  NeutronCandidates::NeutCand GetCurrentLeadingNeutCand() { return fNeutCands.GetMaxCandidate(); };

  NeutronCandidates::NeutCands GetCurrentNeutCands() { return fNeutCands; };

  int GetNNeutCands(){ return fNNeutCands; }
  
 private:
  NeutronCandidates::NeutCands fNeutCands;
  int fNNeutCands;
};

#endif
