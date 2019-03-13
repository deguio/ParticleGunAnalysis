//c++ -o SoverNStudies `root-config --cflags --ldflags --glibs` -lRooFit -lRooFitCore -lRooStats SoverNStudies.cpp SiPM.cc

#include "interface/SiPM.h"
#include "interface/SetTDRStyle.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"


#include <sstream>
#include <algorithm>

#include <TSystem.h>
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


//global objects and variables
TRandom3 *r3      = new TRandom3();

TGraph* outVsIn;
TH2F* biasMap;
int nExp = 0;

/////////////////////////////
/// the actual experiment ///
/////////////////////////////
void runExperiment(CfgManager opts, float nMips, float nPEperMip, float SoN)
{
  //access parameters
  int nEvents = opts.GetOpt<int>("Params.nEvents");

  int nPreciseBins = opts.GetOpt<int>("Params.nPreciseBins");
  float dt = opts.GetOpt<float>("Params.dt");

  float intWin = opts.GetOpt<float>("Params.intWin");

  float tauRise = opts.GetOpt<float>("Params.tauRise");
  float tauDecay = opts.GetOpt<float>("Params.tauDecay");

  int nPixels = opts.GetOpt<int>("Params.nPixels");
  float tau = opts.GetOpt<float>("Params.tauSiPM");
  float xTalk = opts.GetOpt<float>("Params.xTalk");

  float meanN        = std::pow(nPEperMip / SoN, 2); //PE from noise in 25ns
  float meanNtot     = meanN / intWin * dt * nPreciseBins; //PE from noise in the full window

  int minR = 0;
  int maxR = std::max(nMips * nPEperMip, meanN) +100;
  float timeOffset = 0; //100

  //setup SiPM
  HcalSiPM sipm(nPixels, tau, r3);
  sipm.setCrossTalk(xTalk);

  //models and histos
  TF1 sigModel("sigModel","[0]*TMath::Landau(x,[1],[2])", minR, maxR);
  TF1 bkgModel("bkgModel","[0]*TMath::Poisson(x, [1])", minR, maxR);
  TF1 sigTimeModel("sigTimeModel","[0]*(1-exp(-x/[1])) * exp(-x/[2])", 0, nPreciseBins*dt);

  sigModel.SetParameters(1, nMips * nPEperMip, 2);
  sigModel.SetNpx(10000);

  bkgModel.SetParameters(1, meanNtot);  // mean = n^2
  bkgModel.SetNpx(10000);
  bkgModel.SetLineColor(kBlue);


  // TF1* gausSmear = new TF1("gausSmear","[0]*TMath::Gaus(x, [1], [2])", -100, 100);
  // gausSmear->SetParameters(1, 0, gSigma);
  // gausSmear->SetNpx(10000);

  //define time distribution
  sigTimeModel.SetParameters(1, tauRise, tauDecay);
  sigTimeModel.SetNpx(10000);

  //book histos
  int myMax = nMips*nPEperMip*20;
  TH1F h_sig("h_sig","h_sig",myMax,minR,myMax);
  h_sig.GetXaxis()->SetTitle("PE");
  h_sig.SetLineColor(kRed);
  TH1F h_bkg("h_bkg","h_bkg",myMax,minR,myMax);
  h_bkg.GetXaxis()->SetTitle("PE");
  h_bkg.SetLineColor(kBlue);
  TH1F h_tot("h_tot","h_tot",myMax,minR,myMax);
  h_tot.GetXaxis()->SetTitle("PE");
  h_tot.SetLineColor(kBlack);
  h_tot.SetLineWidth(2);

  TProfile arrivalTime("arrivalTime","arrivalTime", nPreciseBins, 0 , nPreciseBins*dt);
  arrivalTime.GetXaxis()->SetTitle("ns");
  TProfile recordedTime("recordedTime","recordedTime", nPreciseBins, 0 , nPreciseBins*dt);
  recordedTime.GetXaxis()->SetTitle("ns");
  recordedTime.SetLineColor(kRed);

  int sumHitsTot = 0;
  int sumPeTot = 0;
  for(int ii=0; ii<nEvents; ++ii)
    {
      //estimate the number of pe associated to a mip traversing the HEback
      float sig = sigModel.GetRandom();
      //estimate the number of pe associated to dark current accordingly to S/N ratio
      float bkg = bkgModel.GetRandom();

      std::vector<unsigned> photonHist(nPreciseBins, 0);
      //account for time distribution signal
      for(int pe = 0; pe < sig; ++pe)
      {
        float t_pe = sigTimeModel.GetRandom();
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


  sumHitsTot /= nEvents;
  sumPeTot /= nEvents;
  outVsIn->SetPoint(nExp, sumPeTot, sumHitsTot);

  ++nExp;
  if(opts.GetOpt<bool>("Params.doMipScan"))
    return;


  //fit smeared signal with langauss
  RooRealVar ee("ee","Signal [pe]",minR,myMax) ;
  //gaussian
  RooRealVar mg("mg","mean gaussian",0) ;
  RooRealVar sg("sg","sigma gaussian",10,0.1,500) ;
  RooGaussian gauss("gauss","gauss",ee,mg,sg) ;
  //landau
  RooRealVar ml("ml","mean landau",nMips*nPEperMip,1,100) ;
  RooRealVar sl("sl","sigma landau",2, 0.5, 500) ;
  RooLandau landau("lx","lx",ee,ml,sl) ;

  ee.setBins(10000,"cache");
  RooFFTConvPdf convsig("convsig","landau (X) gauss",ee,landau,gauss) ;

  RooDataHist dh("dh","dh",ee,Import(h_tot)) ;
  RooPlot* frame = ee.frame(Title("Imported TH1 with Poisson error bars")) ;
  dh.plotOn(frame) ;

  convsig.fitTo(dh,Range(minR,maxR));
  convsig.plotOn(frame) ;

  biasMap->Fill(nPEperMip, SoN, (ml.getVal() - sigModel.GetParameter(1))/sigModel.GetParameter(1));

  TCanvas cPlots;
  cPlots.Divide(2,2);

  cPlots.cd(1);
  h_bkg.GetYaxis()->SetRangeUser(0, std::max(h_bkg.GetMaximum(), h_sig.GetMaximum()));
  h_bkg.DrawCopy();
  h_sig.DrawCopy("sames");
  h_tot.DrawCopy("sames");
  cPlots.Update();

  std::stringstream ss;
  ss << std::fixed << std::setprecision(1) << "S/N: " << SoN << " Signal: " << nPEperMip;
  std::string pars = ss.str();
  TLatex wPoint;
  wPoint.DrawLatexNDC(0.25,0.7, pars.c_str());

  cPlots.cd(2);
  sigModel.GetYaxis()->SetRangeUser(0, std::max(bkgModel.GetMaximum(), sigModel.GetMaximum()));
  sigModel.DrawCopy();
  bkgModel.DrawCopy("sames");
  cPlots.Update();

  cPlots.cd(3);
  frame->Draw();
  cPlots.Update();

  cPlots.cd(4);
  recordedTime.DrawCopy("HISTO");
  arrivalTime.DrawCopy("HISTO,sames");
  cPlots.Update();

  std::stringstream title;
  title << std::fixed << std::setprecision(1) << "plots/nMips_" << nMips << "_SoN_" << SoN << "_Signal_" << nPEperMip << ".pdf";
  std::string oName = title.str();
  cPlots.Print(oName.c_str(),"pdf");
}


////////////
/// main ///
////////////
int main(int argc, char** argv)
{
  gSystem -> Load("CfgManager/lib/libCFGMan.so");
  setTDRStyle();
  gStyle->SetPaintTextFormat("4.3f");
  gROOT->SetBatch(kTRUE);

  CfgManager opts;
  opts.ParseConfigFile(argv[1]);


  float nMipStep = 1.;
  float nMipMin = opts.GetOpt<float>("Params.nMips");
  float nMipMax = nMipMin+nMipStep;
  if(opts.GetOpt<bool>("Params.doMipScan"))
  {
    nMipMin  = opts.GetOpt<float>("MipScan.nMipMin");
    nMipMax  = opts.GetOpt<float>("MipScan.nMipMax");
    nMipStep = opts.GetOpt<float>("MipScan.nMipStep");
  }

  float nPeStep = 1.;
  float nPeMin = opts.GetOpt<float>("Params.nPEperMip");
  float nPeMax = nPeMin+nPeStep;
  if(opts.GetOpt<bool>("Params.doPeScan"))
  {
    nPeMin  = opts.GetOpt<float>("PeScan.nPeMin");
    nPeMax  = opts.GetOpt<float>("PeScan.nPeMax");
    nPeStep = opts.GetOpt<float>("PeScan.nPeStep");
  }

  float SoNStep = 1.;
  float SoNMin = opts.GetOpt<float>("Params.SoN");
  float SoNMax = SoNMin+SoNStep;
  if(opts.GetOpt<bool>("Params.doSoNScan"))
  {
    SoNMin  = opts.GetOpt<float>("SoNScan.SoNMin");
    SoNMax  = opts.GetOpt<float>("SoNScan.SoNMax");
    SoNStep = opts.GetOpt<float>("SoNScan.SoNStep");
  }


  outVsIn   = new TGraph();
  biasMap     = new TH2F("biasMap","biasMap", (nPeMax-nPeMin)/nPeStep, nPeMin-nPeStep/2,nPeMax-nPeStep/2,
                                              (SoNMax-SoNMin)/SoNStep, SoNMin-SoNStep/2,SoNMax-SoNStep/2);

  int count = 1;
  for(float nMips=nMipMin ; nMips < nMipMax; nMips += nMipStep)
    for(float SoN=SoNMin ; SoN < SoNMax; SoN += SoNStep)
      for(float nPEperMip=nPeMin ; nPEperMip < nPeMax; nPEperMip += nPeStep)
      {
        std::cout << count << " --> Running Experiment - MIPS: " << nMips << " SoN: " << SoN << " PEperMIP: " << nPEperMip << std::endl;
        runExperiment(opts, nMips, nPEperMip, SoN);
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
  nonLinFunc->SetParameter(0, opts.GetOpt<int>("Params.nPixels"));
  nonLinFunc->SetNpx(10000);
  nonLinFunc->Draw("sames");

  TF1* nonLinFunc_corr = new TF1("nonLinFunc_corr","[0] * (1 - exp(-x/[0])) / (1 - [1]*exp(-x/[0]))", 0, 50000);
  nonLinFunc_corr->SetParameters(opts.GetOpt<int>("Params.nPixels"), opts.GetOpt<int>("Params.xTalk"));
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


  return 0;
}
