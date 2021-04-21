//File: CVUniverse.h
//Info: Central Value Universe Class. Currently under development.
//      Utilized in following New Systematics Framework (MAT) approach.
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

#include "PlotUtils/MinervaUniverse.h"

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
  virtual int GetNEMBlobs() const { return GetInt("nonvtx_iso_blobs_start_position_z_in_prong_sz"); };
  virtual std::vector<int> GetEMBlobNHitsVec() const { return GetVec<int>("nonvtx_iso_blobs_n_hits_in_prong"); };
  virtual std::vector<double> GetEMBlobEnergyVec() const { return GetVec<double>("nonvtx_iso_blobs_energy_in_prong"); };
  virtual double GetTotalEMBlobNHits() const { 
    std::vector<int> NHitsVec = GetEMBlobNHitsVec();
    return std::accumulate(NHitsVec.begin(), NHitsVec.end(), 0.0);
  };
  virtual double GetTotalEMBlobEnergy(double shift = 0.0) const { 
    std::vector<double> EVec = GetEMBlobEnergyVec();
    return std::accumulate(EVec.begin(), EVec.end(), shift);
  };

  virtual std::vector<double> GetVtx() const { return GetVec<double>("vtx"); };

  virtual int GetNImprovedMichel() const { return GetInt("improved_michel_vertex_type_sz"); };

  virtual int GetNDeadDiscriminatorsUpstreamMuon() const { return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj"); };
  virtual int GetIsMinosMatch() const { return GetInt("isMinosMatchTrack"); };
  virtual int GetNuHelicity() const { return GetInt("MasterAnaDev_nuHelicity"); };

};

#endif
