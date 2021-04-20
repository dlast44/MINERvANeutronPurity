#include "CVUniverse.h"

class NTracksShiftUniverse : public CVUniverse
{
 public:
 NTracksShiftUniverse(typename PlotUtils::MinervaUniverse::config_t chw, double nsigma) : CVUniverse(chw, nsigma)
    {}

  virtual int GetNTracks() const override{
    int shift = 1;
    int shift_val = m_nsigma*shift;
    return shift_val+CVUniverse::GetNTracks();
  }

  virtual std::string ShortName() const { return "NTrackShift"; }
  virtual std::string LongName() const { return "Track No. Shift by 1*nsigma"; }
};
