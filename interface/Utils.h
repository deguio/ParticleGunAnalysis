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
  std::vector<int>* HGCDigiIEta;
  std::vector<int>* HGCDigiIPhi;
  std::vector<int>* HGCDigiLayer;
  std::vector<int>* HGCDigiIndex;

  std::vector<double>* HGCSimHitsIntEnergy;
  std::vector<int>* HGCSimHitsIntIEta;
  std::vector<int>* HGCSimHitsIntIPhi;
  std::vector<int>* HGCSimHitsIntLayer;
  std::vector<int>* HGCSimHitsIntIndex;
  std::vector<double>* HGCSimHitsEta;
  std::vector<double>* HGCSimHitsPhi;

  std::vector<std::vector<float>>* HGCDigiSamples;
};

void InitTreeVars(TTree* chain, MyTreeVars& tt);
TChain* LoadChain(std::string inFileName);


template <class T1, class T2, class T3, class T4>
T1 deltaR2 (T1 eta1, T2 phi1, T3 eta2, T4 phi2) {
  T1 deta = eta1 - eta2;
  T1 dphi = std::abs(phi1-phi2); if (dphi>T1(M_PI)) dphi-=T1(2*M_PI);
  return deta*deta + dphi*dphi;
}

//double undoNonLinearity(double nPix)
//{
//
//  return nPix;
//}

void findPeak(TH1F* h,
  std::vector<std::pair<float,float>>* maxima,
  std::vector<float>* charge,
  std::vector<int> ranges,
  float thr);


void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction);

#endif
