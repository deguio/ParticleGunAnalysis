#ifndef Utils_h
#define Utils_h

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TApplication.h"

#include <utility>
#include <string>
#include <iostream>
#include <cstdlib>


struct MyTreeVars
{
  std::vector<double>* GenParEta;
  std::vector<double>* GenParPhi;
  std::vector<double>* GenParPt;
  std::vector<double>* GenParP;

  std::vector<double>* HGCDigiEta;
  std::vector<double>* HGCDigiPhi;
  std::vector<double>* HGCDigiIEta;
  std::vector<double>* HGCDigiIPhi;
  std::vector<double>* HGCDigiLayer;
  std::vector<double>* HGCDigiIndex;

  std::vector<double>* HGCSimHitsIntEnergy;
  std::vector<double>* HGCSimHitsIntIEta;
  std::vector<double>* HGCSimHitsIntIPhi;
  std::vector<double>* HGCSimHitsIntLayer;
  std::vector<double>* HGCSimHitsIntIndex;
};

void InitTreeVars(TTree* chain, MyTreeVars& tt);
TChain* LoadChain(std::string inFileName);

void findPeak(TH1F* h,
  std::vector<std::pair<float,float>>* maxima,
  std::vector<float>* charge,
  std::vector<int> ranges,
  float thr);


  void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction);

#endif
