#include "interface/Utils.h"


void InitTreeVars(TTree* chain, MyTreeVars& tt)
{
  tt.GenParEta = new std::vector<double>;
  tt.GenParPhi = new std::vector<double>;
  tt.GenParPt = new std::vector<double>;
  tt.GenParP = new std::vector<double>;
  chain -> SetBranchStatus("GenParEta",1); chain -> SetBranchAddress("GenParEta",&tt.GenParEta);
  chain -> SetBranchStatus("GenParPhi",1); chain -> SetBranchAddress("GenParPhi",&tt.GenParPhi);
  chain -> SetBranchStatus("GenParPt",1); chain -> SetBranchAddress("GenParPt",&tt.GenParPt);
  chain -> SetBranchStatus("GenParP",1); chain -> SetBranchAddress("GenParP",&tt.GenParP);

  tt.HGCDigiEta = new std::vector<float>;
  tt.HGCDigiPhi = new std::vector<float>;
  tt.HGCDigiIEta = new std::vector<int>;
  tt.HGCDigiIPhi = new std::vector<int>;
  tt.HGCDigiLayer = new std::vector<int>;
  tt.HGCDigiIndex = new std::vector<int>;
  chain -> SetBranchStatus("HGCDigiEta",1); chain -> SetBranchAddress("HGCDigiEta",&tt.HGCDigiEta);
  chain -> SetBranchStatus("HGCDigiPhi",1); chain -> SetBranchAddress("HGCDigiPhi",&tt.HGCDigiPhi);
  chain -> SetBranchStatus("HGCDigiIEta",1); chain -> SetBranchAddress("HGCDigiIEta",&tt.HGCDigiIEta);
  chain -> SetBranchStatus("HGCDigiIPhi",1); chain -> SetBranchAddress("HGCDigiIPhi",&tt.HGCDigiIPhi);
  chain -> SetBranchStatus("HGCDigiLayer",1); chain -> SetBranchAddress("HGCDigiLayer",&tt.HGCDigiLayer);
  chain -> SetBranchStatus("HGCDigiIndex",1); chain -> SetBranchAddress("HGCDigiIndex",&tt.HGCDigiIndex);

  tt.HGCSimHitsIntEnergy = new std::vector<float>;
  tt.HGCSimHitsIntIEta = new std::vector<int>;
  tt.HGCSimHitsIntIPhi = new std::vector<int>;
  tt.HGCSimHitsIntLayer = new std::vector<int>;
  tt.HGCSimHitsIntIndex = new std::vector<int>;
  tt.HGCSimHitsEta = new std::vector<float>;
  tt.HGCSimHitsPhi = new std::vector<float>;
  chain -> SetBranchStatus("HGCSimHitsIntEnergy",1); chain -> SetBranchAddress("HGCSimHitsIntEnergy",&tt.HGCSimHitsIntEnergy);
  chain -> SetBranchStatus("HGCSimHitsIntIEta",1); chain -> SetBranchAddress("HGCSimHitsIntIEta",&tt.HGCSimHitsIntIEta);
  chain -> SetBranchStatus("HGCSimHitsIntIPhi",1); chain -> SetBranchAddress("HGCSimHitsIntIPhi",&tt.HGCSimHitsIntIPhi);
  chain -> SetBranchStatus("HGCSimHitsIntLayer",1); chain -> SetBranchAddress("HGCSimHitsIntLayer",&tt.HGCSimHitsIntLayer);
  chain -> SetBranchStatus("HGCSimHitsIntIndex",1); chain -> SetBranchAddress("HGCSimHitsIntIndex",&tt.HGCSimHitsIntIndex);
  chain -> SetBranchStatus("HGCSimHitsEta",1); chain -> SetBranchAddress("HGCSimHitsEta",&tt.HGCSimHitsEta);
  chain -> SetBranchStatus("HGCSimHitsPhi",1); chain -> SetBranchAddress("HGCSimHitsPhi",&tt.HGCSimHitsPhi);

  tt.HGCDigiSamples = new std::vector<std::vector<float>>;
  chain -> SetBranchStatus("HGCDigiSamples",1); chain -> SetBranchAddress("HGCDigiSamples",&tt.HGCDigiSamples);
}

//--- load chain ---
TChain* LoadChain(std::string inFileName)
{
  TChain* chain = new TChain("hgcalTupleTree/tree","hgcalTupleTree/tree");

  chain->Add(inFileName.c_str());
  std::cout << " Read " << chain->GetEntries() << " total events in tree " << chain->GetName() << std::endl;

  return chain;
}

//--- find peak ---
void findPeak(TH1F* h,
  std::vector<std::pair<float,float>>* maxima,
  std::vector<float>* charge,
  std::vector<int> ranges,
  float thr)
  {
    std::vector<int> updatedRanges;
    float minMax = -1;
    for(unsigned int rr=1; rr<ranges.size(); rr+=2)
    {
      //to avoid overlapping intervals
      if(ranges.at(rr) <= ranges.at(rr-1))
      continue;

      h->GetXaxis()->SetRange(ranges.at(rr-1),ranges.at(rr));
      float maxY = h->GetMaximum();
      int maxBinX = h->GetMaximumBin();
      float maxX = h->GetXaxis()->GetBinCenter(maxBinX);
      std::pair<float,float> maximum = std::make_pair(maxX,maxY);

      //find the maximum max in this iteration
      if(maxY > minMax)
      minMax = maxY;
      //store only interesting maxima
      if(maxY > thr)
      {
        maxima->push_back(maximum);
        charge->push_back(h->Integral(maxBinX - 12, maxBinX + 12));
      }

      updatedRanges.push_back(ranges.at(rr-1));
      updatedRanges.push_back(maxBinX - 12);
      updatedRanges.push_back(maxBinX + 12);
      updatedRanges.push_back(ranges.at(rr));

    }

    if(minMax < thr)
    return;

    findPeak(h, maxima, charge, updatedRanges, thr);
  }


  /*** find effective sigma ***/
  void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction)
  {
    float integralMax = fraction * histo->Integral();

    int N = histo -> GetNbinsX();
    int M1 = 0;
    int M2 = 0;
    for(int bin1 = 0; bin1 < N; ++bin1)
    {
      if( histo->GetBinContent(bin1+1) > 0. && M1 == 0 ) M1 = bin1-1;
      if( histo->GetBinContent(bin1+1) > 0. ) M2 = bin1+2;
    }

    std::map<int,float> binCenters;
    std::map<int,float> binContents;
    std::map<int,float> binIntegrals;
    for(int bin1 = M1; bin1 < M2; ++bin1)
    {
      binCenters[bin1] = histo->GetBinCenter(bin1+1);
      binContents[bin1] = histo->GetBinContent(bin1+1);

      for(int bin2 = M1; bin2 <= bin1; ++bin2)
      binIntegrals[bin1] += binContents[bin2];
    }

    float min = 0.;
    float max = 0.;
    float delta = 999999.;
    for(int bin1 = M1; bin1 < M2; ++bin1)
    {
      for(int bin2 = bin1+1; bin2 < M2; ++bin2)
      {
        if( (binIntegrals[bin2]-binIntegrals[bin1]) < integralMax ) continue;

        float tmpMin = histo -> GetBinCenter(bin1+1);
        float tmpMax = histo -> GetBinCenter(bin2+1);

        if( (tmpMax-tmpMin) < delta )
        {
          delta = (tmpMax - tmpMin);
          min = tmpMin;
          max = tmpMax;
        }

        break;
      }
    }

    TH1F* smallHisto = (TH1F*)( histo->Clone("smallHisto") );
    for(int bin = 1; bin <= smallHisto->GetNbinsX(); ++bin)
    {
      if( smallHisto->GetBinCenter(bin) < min )
      smallHisto -> SetBinContent(bin,0);

      if( smallHisto->GetBinCenter(bin) > max )
      smallHisto -> SetBinContent(bin,0);
    }
    smallHisto -> SetFillColor(kYellow);
    //smallHisto->Write();

    float mean = smallHisto -> GetMean();
    float meanErr = smallHisto -> GetMeanError();

    ret[0] = mean;
    ret[1] = meanErr;
    ret[2] = min;
    ret[3] = max;
  }

  /*** expo1 ***/
  double expon1(double* x, double* par)
  {
    double xx = x[0];
    double tau = par[0];

    if(xx < 0)
      return 0;
    else
      return (1 - exp(-xx/tau));
  }

  /*** expo2 ***/
  double expon2(double* x, double* par)
  {
    double xx = x[0];
    double tau = par[0];

    return (exp(-xx/tau));
  }


  /*** double expo convolution ***/
  double expoConv(double* x, double* par)
  {
    //[0] = expo1 tau
    //[1] = expo2 tau

    double xx = x[0];
    double tau1 = par[0];
    double tau2 = par[1];

    TF1* expo1 = new TF1("expo1", expon1, -25.,25., 1);
    expo1 -> FixParameter(0,tau1);

    TF1* expo2 = new TF1("expo2", expon2, -25.,25., 1);
    expo2 -> FixParameter(0,tau2);


    // convolute
    double xMin = 0.;
    double xMax = 25.;
    int nSteps = 250;
    double stepWidth = (xMax-xMin)/nSteps;

    double val = 0.;
    for(int i = 0; i < nSteps; ++i)
    {
      double yy = xMin + i*stepWidth;
      val += expo2->Eval(xx) * expo1->Eval(xx-yy);
      //std::cout << "xx = " << xx << " yy = " << yy << " xx-yy = " << xx-yy << std::endl;
    }
    //std::cout << "======" << std::endl;

    delete expo1;
    delete expo2;

    return val;
  }
