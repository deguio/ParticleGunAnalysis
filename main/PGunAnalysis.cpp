#include "interface/SetTDRStyle.h"
#include "interface/Utils.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TApplication.h"


TGraphErrors* makeResGraph(TH2F* resol2d)
{
  TGraphErrors* outGraph = new TGraphErrors();
  int count=0;
  for(int ii = 1; ii < resol2d->GetXaxis()->GetNbins()+1; ++ii)
    {
      //std::cout << "Bin " << ii << " of " << resol2d->GetXaxis()->GetNbins() << std::endl;
      TH1F* tmpHisto = (TH1F*)resol2d->ProjectionY("tmpHisto",ii,ii);
      if(tmpHisto->GetEntries() < 1)
        continue;

      //tmpHisto->Rebin(rebin);
      TF1* mygaus = new TF1("mygauss","gaus",-2, 2);
      tmpHisto->Fit(mygaus,"QR");

      float rangeMin = mygaus->GetParameter(1) - 2*mygaus->GetParameter(2);
      float rangeMax = mygaus->GetParameter(1) + 2*mygaus->GetParameter(2);
      mygaus->SetRange(rangeMin,rangeMax);
      tmpHisto->Fit(mygaus,"QR+");

      float sigma = mygaus->GetParameter(2);
      outGraph->SetPoint(count, resol2d->GetXaxis()->GetBinCenter(ii), sigma);
      outGraph->SetPointError(count, 0, mygaus->GetParError(2));

      // float parstime[4];
      // FindSmallestInterval(parstime, tmpHisto, 0.68);
      // g_timeRes_vsAmp->SetPoint(ii-1,h_timeDiff_vsAmp->GetXaxis()->GetBinCenter(ii), (parstime[3]-parstime[2])/2);

      tmpHisto->Write();
      ++count;

      delete mygaus;
      delete tmpHisto;
    }

  return outGraph;
}


int main(int argc, char** argv)
{
  gSystem -> Load("CfgManager/lib/libCFGMan.so");
  setTDRStyle();

  CfgManager opts;
  opts.ParseConfigFile(argv[1]);

  //ACCESS PARAMETERS
  std::string inputName = opts.GetOpt<std::string>("Input.inputFile");
  std::string outputPath = opts.GetOpt<std::string>("Input.outputFolder");
  std::string outputName = opts.GetOpt<std::string>("Input.outputName");

  int detIndex = opts.GetOpt<int>("Params.detectorIndex");

  float gevPerMip = opts.GetOpt<float>("Params.gevPerMip");
  float mipRange = opts.GetOpt<float>("Params.mipRange");
  int nAdcBits = opts.GetOpt<int>("Params.nAdcBits");
  float zsThr = opts.GetOpt<float>("Params.zsThr");
  int nAdcChannels = pow(2, nAdcBits);

  TApplication *theApp = new TApplication( "app", &argc, argv );

  std::string createFolder = "mkdir -p "+outputPath;
  system(createFolder.c_str());
  std::string outputFile = outputPath+"/"+outputName;
  TFile* outFile = new TFile(outputFile.c_str(), "RECREATE");
  outFile -> cd();

  TChain* chain = LoadChain(inputName);
  MyTreeVars tt;
  InitTreeVars(chain,tt);

  //book histograms
  TH1F* h_simHit_particle_deltaR = new TH1F("h_simHit_particle_deltaR","h_simHit_particle_deltaR",100,0,1);
  TH1F* h_simHit_digi_diff = new TH1F("h_simHit_digi_diff","h_simHit_digi_diff",200,-5,5);
  TH1F* h_simHit_energy = new TH1F("h_simHit_energy","h_simHit_energy",500, 0, 0.005);

  TH1F* h_digi_energy = new TH1F("h_digi_energy","h_digi_energy",150, 0, 10);
  TH1F* h_ped_energy = new TH1F("h_ped_energy","h_ped_energy",150, 0, 10);

  //resolution plots
  TH2F* h_enDiff_vs_ieta = new TH2F("h_enDiff_vs_ieta","h_enDiff_vs_ieta",50, 0, 50, 100, -5, 5);
  TH2F* h_enDiff_vs_layer = new TH2F("h_enDiff_vs_layer","h_enDiff_vs_layer",30, 0, 30, 100, -5, 5);

  //loop over entries
  int nEntries = chain -> GetEntries();
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if(entry % 100 == 0)
    std::cout << "reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    chain -> GetEntry(entry);


    //gen particle
    double genEta = tt.GenParEta->at(0);
    double genPhi = tt.GenParPhi->at(0);

    //loop over ALL SimHits
    for(unsigned int simHit=0; simHit<tt.HGCSimHitsEta->size(); ++simHit)
    {
      double simHitEta = tt.HGCSimHitsEta->at(simHit);
      double simHitPhi = tt.HGCSimHitsPhi->at(simHit);

      double simHit_particle_deltaR = sqrt( deltaR2(genEta,genPhi,simHitEta,simHitPhi) );
      h_simHit_particle_deltaR->Fill(simHit_particle_deltaR);
    }

    //loop over summed simHits
    for(unsigned int simHit=0; simHit<tt.HGCSimHitsIntIEta->size(); ++simHit)
    {
      //select the detector index
      double simHitIndex = tt.HGCSimHitsIntIndex->at(simHit);
      if(simHitIndex != detIndex)
        continue;

      float simHitEn = tt.HGCSimHitsIntEnergy->at(simHit);
      int simHitIEta = tt.HGCSimHitsIntIEta->at(simHit);
      int simHitIPhi = tt.HGCSimHitsIntIPhi->at(simHit);
      int simHitLayer = tt.HGCSimHitsIntLayer->at(simHit);

      //to make a fair comparison with the DIGIs
      if(simHitEn/gevPerMip < zsThr)
        continue;

      h_simHit_energy->Fill(simHitEn);

      for(unsigned int digi=0; digi<tt.HGCDigiIEta->size(); ++digi)
      {
        int digiIEta = tt.HGCDigiIEta->at(digi);
        int digiIPhi = tt.HGCDigiIPhi->at(digi);
        int digiLayer = tt.HGCDigiLayer->at(digi);
        int digiIndex = tt.HGCDigiIndex->at(digi);

        //DetID matching
        if(simHitIEta==digiIEta && simHitIPhi==digiIPhi && simHitLayer==digiLayer && simHitIndex==digiIndex)
        {
          float ped = (tt.HGCDigiSamples->at(digi)[0] + tt.HGCDigiSamples->at(digi)[1])/2;
          float sig = tt.HGCDigiSamples->at(digi)[2] - ped;

          //transform to Mip
          sig = sig / nAdcChannels * mipRange;
          simHitEn = simHitEn / gevPerMip;

          h_digi_energy->Fill(sig);
          h_ped_energy->Fill(ped);

          float diff = sig - simHitEn;
          h_simHit_digi_diff->Fill(diff);
          h_enDiff_vs_ieta->Fill(std::abs(digiIEta), diff);
          h_enDiff_vs_layer->Fill(digiLayer, diff);
          continue;
        }
      }
    }
  }
  std::cout << std::endl;

  //compute resolution
  TDirectory* dirIEta = outFile->mkdir("resVsIEta");
  dirIEta->cd();
  TGraphErrors* g_res_vs_ieta = makeResGraph(h_enDiff_vs_ieta);
  g_res_vs_ieta->SetName("g_res_vs_ieta");
  g_res_vs_ieta->GetXaxis()->SetTitle("IEta");
  g_res_vs_ieta->GetYaxis()->SetTitle("#sigma(reco - gen) [MIP]");

  TCanvas* c_res_vs_ieta = new TCanvas();
  g_res_vs_ieta->Draw("AP");
  c_res_vs_ieta->Update();
  outFile->cd();
  g_res_vs_ieta->Write();

  TDirectory* dirLayer = outFile->mkdir("resVsLayer");
  dirLayer->cd();
  TGraphErrors* g_res_vs_layer = makeResGraph(h_enDiff_vs_layer);
  g_res_vs_layer->SetName("g_res_vs_layer");
  g_res_vs_layer->GetXaxis()->SetTitle("layer");
  g_res_vs_layer->GetYaxis()->SetTitle("#sigma(reco - gen) [MIP]");

  TCanvas* c_res_vs_layer = new TCanvas();
  g_res_vs_layer->Draw("AP");
  c_res_vs_layer->Update();
  outFile->cd();
  g_res_vs_layer->Write();


  //plot
  TCanvas* c_simHit_energy = new TCanvas();
  h_simHit_energy->Draw();
  c_simHit_energy->Update();

  TCanvas* c_simHit_particle_deltaR = new TCanvas();
  h_simHit_particle_deltaR->Draw();
  c_simHit_particle_deltaR->Update();

  TCanvas* c_enDiff_vs_ieta = new TCanvas();
  h_enDiff_vs_ieta->Draw("COLZ");
  c_enDiff_vs_ieta->Update();

  TCanvas* c_enDiff_vs_layer = new TCanvas();
  h_enDiff_vs_layer->Draw("COLZ");
  c_enDiff_vs_layer->Update();


  TCanvas* c_simHit_digi_diff = new TCanvas();
  h_simHit_digi_diff->Draw();
  c_simHit_digi_diff->Update();


  //theApp -> Run();
  outFile->Write();
  outFile->Close();
  return 0;
}
