//c++ -o SoverNStudies `root-config --cflags --ldflags --glibs` -lRooFit -lRooFitCore -lRooStats SoverNStudies.cpp SiPM.cc

#include "interface/SiPM.h"
#include "interface/SetTDRStyle.h"

#include <sstream>
#include <algorithm>

#include <TRandom3.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TApplication.h>
#include <TLatex.h>

#include "RooGlobalFunc.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooFFTConvPdf.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooDataHist.h"

using namespace RooFit ;

//global parameters
int   nEvents      = 1000;
int   nMips        = 1;

//experiments
float SoN_min = 1;
float SoN_max = 15; //15
float SoN_step = 1;

float nPE_min = 5;
float nPE_max = 35; //35
float nPE_step = 5;


//time granularity
int   nPreciseBins = 250; //1250
float dt           = 0.1; //ns

int   nPEperMip    = 30; //30 nPE per MIP
float SoN          = 10;  //integration window for noise = 25ns
float intWin       = 25;
//float gSigma       = 0.05;
float meanN        = std::pow(nPEperMip / SoN, 2); //PE from noise in 25ns
float meanNtot     = meanN / intWin * dt * nPreciseBins; //PE from noise in the full window
int   maxMips      = 1000; //69 saturate the ADC (just to have 1 MIP around ~15ADC)
int   nBits        = 10;

//scintillation shape
float tauRise      = 0.9;   //0.9
float tauDecay     = 2.1; //2.1

//SiPM params
int nPixels        = 7500 ; //7500
float tau          = 10;   //10
float xTalk        = 0.2;    //0.2

//int nBins = std::pow(2, nBits);
int nBins = maxMips * nPEperMip; //for now just plot the nPE
int minR = 0;
int maxR = maxMips * nPEperMip;
float timeOffset = 0; //100

//global objects and variables
TRandom3 *r3      = new TRandom3();

TF1* sigModel     = new TF1("sigModel","[0]*TMath::Landau(x,[1],[2])", minR, maxR);
TF1* bkgModel     = new TF1("bkgModel","[0]*TMath::Poisson(x, [1])", minR, maxR);
TF1* sigTimeModel = new TF1("sigTimeModel","[0]*(1-exp(-x/[1])) * exp(-x/[2])", 0, nPreciseBins*dt);

TGraph* outVsIn   = new TGraph();
TH2F* biasMap     = new TH2F("biasMap","biasMap", (nPE_max-nPE_min)/nPE_step, nPE_min-nPE_step/2,nPE_max-nPE_step/2,
                                                  (SoN_max-SoN_min)/SoN_step, SoN_min-SoN_step/2,SoN_max-SoN_step/2);
int nExp = 0;

/////////////////////////////
/// the actual experiment ///
/////////////////////////////
void runExperiment()
{
  meanN    = std::pow(nPEperMip / SoN, 2);
  meanNtot = meanN / intWin * dt * nPreciseBins;
  nBins    = maxMips * nPEperMip; //for now just plot the nPE
  maxR     = std::max(maxMips * nPEperMip, (int)(meanN*1.5));

  sigModel->SetRange(minR, maxR);
  bkgModel->SetRange(minR, maxR);

  HcalSiPM sipm(nPixels, tau, r3);
  sipm.setCrossTalk(xTalk);

  sigModel->SetParameters(1, nMips*nPEperMip, 2);
  sigModel->SetNpx(10000);

  bkgModel->SetParameters(1, meanNtot);  // mean = n^2
  bkgModel->SetNpx(10000);
  bkgModel->SetLineColor(kBlue);


  // TF1* gausSmear = new TF1("gausSmear","[0]*TMath::Gaus(x, [1], [2])", -100, 100);
  // gausSmear->SetParameters(1, 0, gSigma);
  // gausSmear->SetNpx(10000);

  //define time distribution
  sigTimeModel->SetParameters(1, tauRise, tauDecay);
  sigTimeModel->SetNpx(10000);

  //book histos
  int myMax = nMips*nPEperMip*20;
  TH1F h_sig("h_sig","h_sig",myMax,minR,myMax);
  h_sig.GetXaxis()->SetTitle("PE");
  h_sig.SetLineColor(kRed);
  TH1F h_bkg("h_bkg","h_bkg",myMax,minR,myMax);
  h_bkg.SetLineColor(kBlue);
  TH1F h_tot("h_tot","h_tot",myMax,minR,myMax);
  h_tot.SetLineColor(kBlack);
  h_tot.SetLineWidth(2);

  TProfile arrivalTime("arrivalTime","arrivalTime", nPreciseBins, 0 , nPreciseBins*dt);
  TProfile recordedTime("recordedTime","recordedTime", nPreciseBins, 0 , nPreciseBins*dt);
  recordedTime.SetLineColor(kRed);

  int sumHitsTot = 0;
  int sumPeTot = 0;
  for(int ii=0; ii<nEvents; ++ii)
    {
      //estimate the number of pe associated to a mip traversing the HEback
      float sig = sigModel->GetRandom();
      //estimate the number of pe associated to dark current accordingly to S/N ratio
      float bkg = bkgModel->GetRandom();

      std::vector<unsigned> photonHist(nPreciseBins, 0);
      //account for time distribution signal
      for(int pe = 0; pe < sig; ++pe)
      {
        float t_pe = sigTimeModel->GetRandom();
        int t_bin = (int)((t_pe + timeOffset) / dt);
        if(t_bin >= 0 and (static_cast<unsigned>(t_bin) < photonHist.size()))
          photonHist[t_bin] += 1;
      }
      //account for time distribution bkg
      for(int pe = 0; pe < bkg; ++pe)
      {
        int t_bin = (int)(nPreciseBins * r3->Rndm());
        if(t_bin >= 0 and (static_cast<unsigned>(t_bin) < photonHist.size()))
          photonHist[t_bin] += 1;
      }
      //account for time distribution background
      // for(unsigned tbin = 0; tbin < photonHist.size(); ++tbin)
      // {
      //   int dcNPE = std::floor(bkgModel->GetRandom()) ;
      //   photonHist[tbin] += dcNPE;
      // }

      // float smear = 0;
      // for (int ii=0; ii<tot; ++ii)
      //   {
	    //     smear += gausSmear->GetRandom();
      //   }

      //reset SiPM
      sipm.setNCells(nPixels);

      //evaluate SiPM response
      double elapsedTime = 0.;
      unsigned sumPE = 0;
      double sumHits = 0;
      for(unsigned tbin = 0; tbin < photonHist.size(); ++tbin)
      {
				unsigned pe = photonHist[tbin];
        double hitted = 0;
				if(pe > 0)
        {
          hitted = sipm.hitCells(r3, pe, 0., elapsedTime);

          if(tbin*dt >= timeOffset && tbin*dt < timeOffset+intWin)
          {
            sumPE += pe;
            sumHits += hitted;

            sumPeTot += pe;
            sumHitsTot += hitted;
          }
				}
        arrivalTime.Fill(elapsedTime, pe);
        recordedTime.Fill(elapsedTime, hitted);

				elapsedTime += dt;
			}

      h_sig.Fill(sig);
      float inWinBkg = bkg / (nPreciseBins * dt / intWin);
      h_bkg.Fill(inWinBkg);
      h_tot.Fill(sig+bkg - meanN);
    }

  //fit smeared signal with langauss
  RooRealVar* ee = new RooRealVar("ee","Signal [pe]",minR,myMax) ;
  //gaussian
  RooRealVar* mg = new RooRealVar("mg","mean gaussian",0) ;
  RooRealVar* sg = new RooRealVar("sg","sigma gaussian",10,0.1,500) ;
  RooAbsPdf* gauss = new RooGaussian("gauss","gauss",*ee,*mg,*sg) ;
  //landau
  RooRealVar* ml = new RooRealVar("ml","mean landau",nMips*nPEperMip,1,100) ;
  RooRealVar* sl = new RooRealVar("sl","sigma landau",2, 0.5, 500) ;
  RooAbsPdf* landau = new RooLandau("lx","lx",*ee,*ml,*sl) ;

  ee->setBins(10000,"cache");
  RooAbsPdf* convsig = new RooFFTConvPdf("convsig","landau (X) gauss",*ee,*landau,*gauss) ;

  RooDataHist dh("dh","dh",*ee,Import(h_tot)) ;
  RooPlot* frame = ee->frame(Title("Imported TH1 with Poisson error bars")) ;
  dh.plotOn(frame) ;

  // convsig->fitTo(dh,Range(minR,maxR));
  // convsig->plotOn(frame) ;

  biasMap->Fill(nPEperMip, SoN, (ml->getVal() - sigModel->GetParameter(1))/sigModel->GetParameter(1));

  sumHitsTot /= nEvents;
  sumPeTot /= nEvents;
  outVsIn->SetPoint(nExp, sumPeTot, sumHitsTot);
  ++nExp;

  TCanvas* cPlots = new TCanvas();
  cPlots->Divide(2,2);

  cPlots->cd(1);
  h_bkg.GetYaxis()->SetRangeUser(0, std::max(h_bkg.GetMaximum(), h_sig.GetMaximum()));
  h_bkg.DrawCopy();
  h_sig.DrawCopy("sames");
  h_tot.DrawCopy("sames");
  cPlots->Update();

  std::stringstream ss;
  ss << std::fixed << std::setprecision(1) << "S/N: " << SoN << " Signal: " << nPEperMip;
  std::string pars = ss.str();
  TLatex wPoint;
  wPoint.DrawLatexNDC(0.25,0.7, pars.c_str());

  cPlots->cd(2);
  bkgModel->GetYaxis()->SetRangeUser(0, std::max(bkgModel->GetMaximum(), sigModel->GetMaximum()));
  bkgModel->Draw();
  sigModel->Draw("sames");
  cPlots->Update();

  cPlots->cd(3);
  frame->Draw();
  cPlots->Update();

  cPlots->cd(4);
  arrivalTime.DrawCopy("HISTO");
  recordedTime.DrawCopy("HISTO,sames");
  cPlots->Update();

  std::stringstream title;
  title << std::fixed << std::setprecision(1) << "plots/SoN_" << SoN << "_Signal_" << nPEperMip << ".pdf";
  std::string oName = title.str();
  cPlots->Print(oName.c_str(),"pdf");
}


////////////
/// main ///
////////////
int main(int argc, char** argv)
{
  TApplication *theApp = new TApplication( "app", &argc, argv );
  setTDRStyle();
  gStyle->SetPaintTextFormat("4.3f");
  gROOT->SetBatch(kTRUE);

  int count = 1;
  for( ; nMips < 1000; nMips += 100)
  // for(SoN=SoN_min ; SoN < SoN_max; SoN += SoN_step)
  //   for(nPEperMip=nPE_min ; nPEperMip < nPE_max; nPEperMip += nPE_step)
    {
      std::cout << count << " --> Running Experiment - MIPS: " << nMips << " SoN: " << SoN << " PEperMIP: " << nPEperMip << std::endl;
      runExperiment();
      ++count;
    }

  TCanvas* cLin = new TCanvas("cLin","cLin");
  cLin->SetGridx();
  cLin->SetGridy();

  outVsIn->SetMarkerStyle(20);
  outVsIn->Draw("AP");

  outVsIn->GetYaxis()->SetTitle("output [pe]");
  outVsIn->GetXaxis()->SetTitle("input [pe]");

  TF1* nonLinFunc = new TF1("nonLinFunc","[0]*(1 - exp(-x/[0]))", 0, 50000);
  nonLinFunc->SetParameter(0, nPixels);
  nonLinFunc->SetNpx(10000);
  nonLinFunc->Draw("sames");

  TF1* nonLinFunc_corr = new TF1("nonLinFunc_corr","[0] * (1 - exp(-x/[0])) / (1 - [1]*exp(-x/[0]))", 0, 50000);
  nonLinFunc_corr->SetParameters(nPixels, xTalk);
  nonLinFunc_corr->SetNpx(10000);
  nonLinFunc_corr->SetLineColor(kBlue);
  nonLinFunc_corr->Draw("sames");

  cLin->Update();
  cLin->Print("plots/lin.pdf","pdf");


  TCanvas* cBias = new TCanvas("cBias","cBias");
  biasMap->GetXaxis()->SetTitle("Signal [pe]");
  biasMap->GetYaxis()->SetTitle("S/N");
  biasMap->GetZaxis()->SetTitle("bias");
  biasMap->Draw("COLZ,TEXT");
  cBias->Update();

  cBias->Print("plots/bias.pdf","pdf");


//  theApp -> Run();
  return 0;
}
