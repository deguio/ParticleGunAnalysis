#include "interface/SiPM.h"
#include "interface/SetTDRStyle.h"
#include "interface/Utils.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"


#include <sstream>
#include <algorithm>

#include <TSystem.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile2D.h>
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
TFile* outFile;

TGraph* outVsIn;
TGraph* inVsOut;
TH2F* biasMap;
TH2F* biasErrorMap;
TProfile2D* resMap;
int nExp = 0;

/////////////////////////////
/// the actual experiment ///
/////////////////////////////
void runExperiment(CfgManager opts, float nMips, float nPEperMip, float SoN)
{
  std::string folder = "nMips_"+std::to_string(nMips)+"_SoN_"+std::to_string(SoN)+"_Signal_"+std::to_string(nPEperMip);
  TDirectory* rootFolder = outFile->mkdir(folder.c_str());
  rootFolder->cd();

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
  int maxR = std::max(nMips * nPEperMip, meanN) + 300;
  float timeOffset = 0; //100

  //setup SiPM
  HcalSiPM sipm(nPixels, tau, r3);
  sipm.setCrossTalk(xTalk);

  //models and histos
  TF1 sigModel("sigModel","[0]*TMath::Landau(x,[1],[2])", minR, maxR);
  TF1 sigSmearModel("sigSmearModel","[0]*TMath::Poisson(x, [1])", minR, maxR);
  TF1 bkgModel("bkgModel","[0]*TMath::Poisson(x, [1])", minR, maxR);
  TF1 sigTimeModel("sigTimeModel", expoConv, 0., nPreciseBins*dt, 2);

  sigModel.SetParameters(1, nMips * nPEperMip, 2);
  sigModel.SetNpx(10000);

  bkgModel.SetParameters(1, meanNtot);  // mean = n^2
  bkgModel.SetNpx(10000);
  bkgModel.SetLineColor(kBlue);


  // TF1* gausSmear = new TF1("gausSmear","[0]*TMath::Gaus(x, [1], [2])", -100, 100);
  // gausSmear->SetParameters(1, 0, gSigma);
  // gausSmear->SetNpx(10000);

  //define time distribution
  sigTimeModel.SetParameters(tauRise, tauDecay);
  //sigTimeModel.SetNpx(10000);

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

  TH1F h_res("h_res","h_res", 4400, -5, 50);

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
      //smear the number of PE from the signal deposit according to poisson
      sigSmearModel.SetParameters(1., sig);
      sig = sigSmearModel.GetRandom();

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
      h_res.Fill((sig+bkg - meanN - nPEperMip) / nPEperMip);
    }


  sumHitsTot /= nEvents;
  sumPeTot /= nEvents;
  outVsIn->SetPoint(nExp, sumPeTot, sumHitsTot);
  inVsOut->SetPoint(nExp, sumHitsTot, (float)sumPeTot/(float)sumHitsTot);

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

  float ret[4];
  FindSmallestInterval(ret, &h_tot, 0.68);
  float range = ret[3] - ret[2];

  convsig.fitTo(dh,Range(minR, range*10));
  convsig.plotOn(frame) ;

  int bin = biasMap->FindBin(nPEperMip, SoN);
  biasMap->SetBinContent(bin, (ml.getVal() - sigModel.GetParameter(1))/sigModel.GetParameter(1));
  biasMap->SetBinError(bin, ml.getError()/sigModel.GetParameter(1));
  biasErrorMap->SetBinContent(bin, ml.getError()/sigModel.GetParameter(1));

  //fit res plots
//  TF1* mygaus = new TF1("mygauss","gaus",-2, 2);
//  h_res.Fit(mygaus,"QR");
//  float rangeMin = mygaus->GetParameter(1) - 2*mygaus->GetParameter(2);
//  float rangeMax = mygaus->GetParameter(1) + 2*mygaus->GetParameter(2);
//  mygaus->SetRange(rangeMin,rangeMax);
//  h_res.Fit(mygaus,"QR+");
//
//  float sigma = mygaus->GetParameter(2);
//  float sigmaE = mygaus->GetParError(2);

  float ret2[4];
  FindSmallestInterval(ret2, &h_res, 0.68);
  float range2 = ret2[3] - ret2[2];

  resMap->Fill(nPEperMip, SoN, range2/2.);


  //plot
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

  TCanvas cRes;
  h_res.DrawCopy();
  cRes.Update();


  std::string outFolder = opts.GetOpt<std::string>("Input.outputFolder");
  std::stringstream title;
  title << std::fixed << std::setprecision(1) << outFolder << "/nMips_" << nMips << "_SoN_" << SoN << "_Signal_" << nPEperMip << ".pdf";
  std::string oName = title.str();
  cPlots.Print(oName.c_str(),"pdf");
  frame->Write("langaus");

  std::stringstream title2;
  title2 << std::fixed << std::setprecision(1) << outFolder << "/hRes_" << nMips << "_SoN_" << SoN << "_Signal_" << nPEperMip << ".pdf";
  oName = title2.str();
  cRes.Print(oName.c_str(),"pdf");
  h_res.Write();

  outFile->cd();

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
  std::string outFolder = opts.GetOpt<std::string>("Input.outputFolder");

  std::string createFolder = "mkdir -p "+outFolder;
  system(createFolder.c_str());


  outFile = new TFile((outFolder+"/plots.root").c_str(), "RECREATE");
  outFile -> cd();


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
  inVsOut   = new TGraph();
  biasMap     = new TH2F("biasMap","biasMap", (nPeMax-nPeMin)/nPeStep, nPeMin-nPeStep/2,nPeMax-nPeStep/2,
                                              (SoNMax-SoNMin)/SoNStep, SoNMin-SoNStep/2,SoNMax-SoNStep/2);
  biasErrorMap     = new TH2F("biasErrorMap","biasErrorMap", (nPeMax-nPeMin)/nPeStep, nPeMin-nPeStep/2,nPeMax-nPeStep/2,
                                                             (SoNMax-SoNMin)/SoNStep, SoNMin-SoNStep/2,SoNMax-SoNStep/2);

  resMap  = new TProfile2D("resMap","resMap", (nPeMax-nPeMin)/nPeStep, nPeMin-nPeStep/2,nPeMax-nPeStep/2,
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

  //LINEARITY
  TCanvas* cLin = new TCanvas("cLin","cLin");
  cLin->SetGridx();
  cLin->SetGridy();

  outVsIn->SetMarkerStyle(20);
  outVsIn->Draw("AP");
  outVsIn->GetYaxis()->SetTitle("output [pe]");
  outVsIn->GetXaxis()->SetTitle("input [pe]");

  TF1* nonLinFunc = new TF1("nonLinFunc","[0]*(1 - exp(-x/[0]))", 0, 80000);
  nonLinFunc->SetParameter(0, opts.GetOpt<int>("Params.nPixels"));
  nonLinFunc->SetNpx(10000);
  nonLinFunc->Draw("sames");
  cLin->Update();

  TF1* nonLinFunc_corr = new TF1("nonLinFunc_corr","[0] * (1 - exp(-x/[0])) / (1 - [1]*exp(-x/[0]))", 0, 80000);
  nonLinFunc_corr->SetParameters(opts.GetOpt<int>("Params.nPixels"), opts.GetOpt<float>("Params.xTalk"));
  nonLinFunc_corr->SetNpx(10000);
  nonLinFunc_corr->SetLineColor(kBlue);
  nonLinFunc_corr->Draw("sames");
  cLin->Update();
  cLin->Print((outFolder+"/lin.pdf").c_str(),"pdf");

  //CORRECTION
  TCanvas* cCorr = new TCanvas("cCorr","cCorr");
  cCorr->SetGridx();
  cCorr->SetGridy();

  inVsOut->SetMarkerStyle(20);
  inVsOut->Draw("AP");
  inVsOut->GetYaxis()->SetTitle("corr");
  inVsOut->GetXaxis()->SetTitle("output [pe]");

  TF1* corrFunc = new TF1("corrFunc","pol2", 0, 50000);
  inVsOut->Fit(corrFunc,"R");
  cCorr->Update();
  cCorr->Print((outFolder+"/corr.pdf").c_str(),"pdf");


  //BIAS
  TCanvas* cBias = new TCanvas("cBias","cBias");
  biasMap->GetXaxis()->SetTitle("Signal [pe]");
  biasMap->GetYaxis()->SetTitle("S/N");
  biasMap->GetZaxis()->SetTitle("bias");
  biasMap->Draw("COLZ,TEXT");
  cBias->Update();
  cBias->Print((outFolder+"/bias.pdf").c_str(),"pdf");

  TCanvas* cBiasError = new TCanvas("cBiasError","cBiasError");
  biasErrorMap->GetXaxis()->SetTitle("Signal [pe]");
  biasErrorMap->GetYaxis()->SetTitle("S/N");
  biasErrorMap->GetZaxis()->SetTitle("relative error on landau location par");
  biasErrorMap->Draw("COLZ,TEXT");
  cBiasError->Update();
  cBiasError->Print((outFolder+"/biasError.pdf").c_str(),"pdf");

  //RES
  TCanvas* cRes = new TCanvas("cRes","cRes");
  resMap->GetXaxis()->SetTitle("Signal [pe]");
  resMap->GetYaxis()->SetTitle("S/N");
  resMap->GetZaxis()->SetTitle("Resolution");
  resMap->Draw("COLZ,TEXT");
  cRes->Update();

  cRes->Print((outFolder+"/resolution.pdf").c_str(),"pdf");

  outFile->Write();
  outFile->Close();
  return 0;
}
