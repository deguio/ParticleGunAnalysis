#include "interface/SetTDRStyle.h"
#include "interface/Utils.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include "TSystem.h"
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



int main(int argc, char** argv)
{
  gSystem -> Load("CfgManager/lib/libCFGMan.so");
  setTDRStyle();

  CfgManager opts;
  opts.ParseConfigFile(argv[1]);

  std::string inputName = opts.GetOpt<std::string>("Input.inputFile");
  std::string outputPath = opts.GetOpt<std::string>("Input.outputFolder");

  TApplication *theApp = new TApplication( "app", &argc, argv );

  std::string createFolder = "mkdir -p "+outputPath;
  system(createFolder.c_str());
  std::string outputFile = outputPath+"/PGunAnalysis.root";
  TFile* outFile = new TFile(outputFile.c_str(), "RECREATE");
  outFile -> cd();

  TChain* chain = LoadChain(inputName);
  MyTreeVars tt;
  InitTreeVars(chain,tt);

  //book histograms
  TH1F* h_simHit_particle_deltaR = new TH1F("h_simHit_particle_deltaR","h_simHit_particle_deltaR",100,0,2);
  TH1F* h_simHit_digi_diff = new TH1F("h_simHit_digi_diff","h_simHit_digi_diff",200,-50,50);


  int nEntries = chain -> GetEntries();
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if(entry % 1000 == 0)
    std::cout << "reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    chain -> GetEntry(entry);

    double genEta = tt.GenParEta->at(0);
    double genPhi = tt.GenParPhi->at(0);

    //loop over all SimHits
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
      //take only scintillators
      double simHitIndex = tt.HGCSimHitsIntIndex->at(simHit);
      if(simHitIndex != 2)
        continue;

      float simHitEn = tt.HGCSimHitsIntEnergy->at(simHit);
      int simHitIEta = tt.HGCSimHitsIntIEta->at(simHit);
      int simHitIPhi = tt.HGCSimHitsIntIPhi->at(simHit);
      int simHitLayer = tt.HGCSimHitsIntLayer->at(simHit);
      for(unsigned int digi=0; digi<tt.HGCDigiIEta->size(); ++digi)
      {
        int digiIEta = tt.HGCDigiIEta->at(digi);
        int digiIPhi = tt.HGCDigiIPhi->at(digi);
        int digiLayer = tt.HGCDigiLayer->at(digi);
        int digiIndex = tt.HGCDigiIndex->at(digi);
        if(simHitIEta==digiIEta && simHitIPhi==digiIPhi && simHitLayer==digiLayer && simHitIndex==digiIndex)
        {
          float ped = (tt.HGCDigiSamples->at(digi)[1] + tt.HGCDigiSamples->at(digi)[2])/2;
          float sig = tt.HGCDigiSamples->at(digi)[9] - ped;

          // 1 MIP = 500 KeV = 0.0005 GeV
          double res = (sig - simHitEn/0.0005)/(simHitEn/0.0005);

          h_simHit_digi_diff->Fill(res);
          continue;
        }
      }
    }



  }
  std::cout << std::endl;


  TCanvas* c_simHit_particle_deltaR = new TCanvas();
  h_simHit_particle_deltaR->Draw();
  c_simHit_particle_deltaR->Update();

  TCanvas* c_simHit_digi_diff = new TCanvas();
  h_simHit_digi_diff->Draw();
  c_simHit_digi_diff->Update();

  theApp -> Run();
  outFile->Write();
  outFile->Close();
  return 0;
}
