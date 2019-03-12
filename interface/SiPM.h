/**
  \class HcalSiPM
  \brief A general implementation for the response of a SiPM.
*/
#include <vector>
#include <algorithm>
#include <unordered_map>

#include <TRandom3.h>

class HcalSiPM {
 public:
  HcalSiPM(int nCells = 1, double tau = 15., TRandom3* rnd = 0);

  virtual ~HcalSiPM();

  void resetSiPM() { std::fill(theSiPM.begin(), theSiPM.end(), -999.); }
  virtual double hitCells(TRandom3* , unsigned int pes, double tempDiff = 0., double photonTime = 0.);


  virtual double totalCharge() const { return totalCharge(theLastHitTime); }
  virtual double totalCharge(double time) const;

  int    getNCells()      const { return theCellCount; }
  double getTau()         const { return theTau; }
  double getCrossTalk()   const { return theCrossTalk; }
  double getTempDep()     const { return theTempDep; }

  void setNCells(int nCells);
  void setTau(double tau);
  void setCrossTalk(double xtalk); //  Borel-Tanner "lambda"
  void setTemperatureDependence(double tempDep);
  void setSaturationPars(const std::vector<float>& pars);

 protected:

  typedef std::pair<unsigned int, std::vector<double> > cdfpair;
  typedef std::unordered_map< unsigned int, cdfpair > cdfmap;

  // void expRecover(double dt);

  double cellCharge(double deltaTime) const;
  unsigned int addCrossTalkCells(TRandom3* engine, unsigned int in_pes);

  //numerical random generation from Borel-Tanner distribution
  double Borel(unsigned int n, double lambda, unsigned int k);
  const cdfpair& BorelCDF(unsigned int k);

  unsigned int theCellCount;
  std::vector< double > theSiPM;
  double theTau;
  double theTauInv;
  double theCrossTalk;
  double theTempDep;
  double theLastHitTime;

  cdfmap borelcdfs;

  TRandom3* engine;
};
