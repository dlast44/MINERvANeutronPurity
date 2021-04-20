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
  virtual int GetNNeutBlobs() const { return GetInt("MasterAnaDev_BlobIs3D_sz");};

};

#endif
