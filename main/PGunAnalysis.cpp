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

//                                             minVal,   maxVal,   nBinsToMerge
TGraphErrors* makeResGraph(TH2F* resol2d, float min, float max, int width)
{
  TGraphErrors* outGraph = new TGraphErrors();
  int count=0;
  for(int ii = resol2d->GetXaxis()->FindBin(min); ii < resol2d->GetXaxis()->FindBin(max)+1; ii+=width)
  {
    //std::cout << "Bin " << ii << " of " << resol2d->GetXaxis()->GetNbins() << std::endl;
    TH1F* tmpHisto = (TH1F*)resol2d->ProjectionY((std::string(std::to_string(ii)+"_to_"+std::to_string(ii+width))).c_str(),ii,ii + width-1);
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
    float sigmaE = mygaus->GetParError(2);
    outGraph->SetPoint(count, resol2d->GetXaxis()->GetBinCenter(ii), sigma);
    outGraph->SetPointError(count, 0, sigmaE);

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
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  CfgManager opts;
  opts.ParseConfigFile(argv[1]);

  int nEntries = -1;
  if(argc > 2)
  nEntries = atoi(argv[2]);

  //ACCESS PARAMETERS
  std::string inputName = opts.GetOpt<std::string>("Input.inputFile");
  std::string outputPath = opts.GetOpt<std::string>("Input.outputFolder");
  std::string outputName = opts.GetOpt<std::string>("Input.outputName");

  int detIndex = opts.GetOpt<int>("Params.detectorIndex");

  float gevPerMip = opts.GetOpt<float>("Params.gevPerMip");
  float mipRange = opts.GetOpt<float>("Params.mipRange");
  int nAdcBits = opts.GetOpt<int>("Params.nAdcBits");
  float zsThr = opts.GetOpt<float>("Params.zsThr");
  std::vector<double> radiusRegions = opts.GetOpt<std::vector<double>>("Plot.radiusRegions");
  int nAdcChannels = pow(2, nAdcBits);

  //TApplication *theApp = new TApplication( "app", &argc, argv );

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
  TH1F* h_simHit_energy = new TH1F("h_simHit_energy","h_simHit_energy",150, 0, 10);

  TH1F* h_digi_energy = new TH1F("h_digi_energy","h_digi_energy",150, 0, 10);
  TH1F* h_ped_energy = new TH1F("h_ped_energy","h_ped_energy",150, 0, 10);

  TProfile* h_pulseShape = new TProfile("h_pulseShape","h_pulseShape", 5, 0., 5.);

  //resolution plots
  TH2F* h_enDiff_vs_irad = new TH2F("h_enDiff_vs_irad","h_enDiff_vs_irad",50, 0, 50, 200, -5, 5);
  TH2F* h_enDiff_vs_layer = new TH2F("h_enDiff_vs_layer","h_enDiff_vs_layer",30, 0, 30, 100, -5, 5);
  std::map<double, TH2F*> hMap_enDiff_vs_layer;
  for(unsigned int bin=0; bin<radiusRegions.size()-1; ++bin)
    hMap_enDiff_vs_layer[radiusRegions[bin]] = new TH2F(Form("h_enDiff_vs_layer_%f-%f",radiusRegions[bin],radiusRegions[bin+1]),"",30, 0, 30, 100, -5, 5);

  //loop over entries
  if(nEntries == -1)
  nEntries = chain -> GetEntries();
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if(entry % 100 == 0)
    std::cout << "reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    chain -> GetEntry(entry);


    //gen particle
    double genEta = tt.GenParEta->at(0);
    double genPhi = tt.GenParPhi->at(0);

    //loop over ALL SimHits
    //    for(unsigned int simHit=0; simHit<tt.HGCSimHitsEta->size(); ++simHit)
    //    {
    //      double simHitEta = tt.HGCSimHitsEta->at(simHit);
    //      double simHitPhi = tt.HGCSimHitsPhi->at(simHit);
    //
    //      double simHit_particle_deltaR = sqrt( deltaR2(genEta,genPhi,simHitEta,simHitPhi) );
    //      h_simHit_particle_deltaR->Fill(simHit_particle_deltaR);
    //    }

    //loop over digis
    for(unsigned int digi=0; digi<tt.HGCDigiIEta->size(); ++digi)
    {
      int digiIndex = tt.HGCDigiIndex->at(digi);
      //consider only the relevant subdet
      if(digiIndex != detIndex)
        continue;

      int digiIRad = tt.HGCDigiIEta->at(digi);
      int digiIPhi = tt.HGCDigiIPhi->at(digi);
      int digiLayer = tt.HGCDigiLayer->at(digi);
      float ped = (tt.HGCDigiSamples->at(digi)[0] + tt.HGCDigiSamples->at(digi)[1])/2;
      float sig = tt.HGCDigiSamples->at(digi)[2];

      //transform to Mips
      sig = sig / nAdcChannels * mipRange;
      ped = ped / nAdcChannels * mipRange;

      //apply ZS on digis before subtracting ped
      if(sig <= zsThr)
        continue;

      //subtract the ped
      sig = sig - ped;

      //loop over summed simHits
      for(unsigned int simHit=0; simHit<tt.HGCSimHitsIntIEta->size(); ++simHit)
      {
        //select the detector index
        double simHitIndex = tt.HGCSimHitsIntIndex->at(simHit);
        if(simHitIndex != detIndex)
          continue;

        int simHitIRad = tt.HGCSimHitsIntIEta->at(simHit);
        int simHitIPhi = tt.HGCSimHitsIntIPhi->at(simHit);
        int simHitLayer = tt.HGCSimHitsIntLayer->at(simHit);
        float simHitEn = tt.HGCSimHitsIntEnergy->at(simHit);

        //transform to Mips
        simHitEn = simHitEn / gevPerMip;

        //DetID matching
        if(simHitIRad==digiIRad && simHitIPhi==digiIPhi && simHitLayer==digiLayer)
        {
          //fill pulse shape plot
          for(unsigned int sampl = 0; sampl<tt.HGCDigiSamples->at(digi).size(); ++sampl)
            h_pulseShape->Fill(sampl, tt.HGCDigiSamples->at(digi)[sampl]);


          h_digi_energy->Fill(sig);
          h_ped_energy->Fill(ped);
          h_simHit_energy->Fill(simHitEn);

          float diff = sig - simHitEn;
          h_simHit_digi_diff->Fill(diff);
          h_enDiff_vs_irad->Fill(std::abs(digiIRad), diff/simHitEn);
          h_enDiff_vs_layer->Fill(digiLayer,         diff/simHitEn);

          //fill histo map
          for(auto itr=hMap_enDiff_vs_layer.rbegin(); itr!=hMap_enDiff_vs_layer.rend(); ++itr)
          {
            if(std::abs(digiIRad)>itr->first)
            {
              itr->second->Fill(digiLayer, diff/simHitEn);
              break;
            }
          }

          break;
        }
      }
    }
  }
  std::cout << std::endl;

  //compute resolution
  TDirectory* dirIRad = outFile->mkdir("resVsIRad");
  dirIRad->cd();
  TGraphErrors* g_res_vs_irad = makeResGraph(h_enDiff_vs_irad, 0, 50, 2);
  g_res_vs_irad->SetName("g_res_vs_irad");
  g_res_vs_irad->GetXaxis()->SetTitle("IRad");
  g_res_vs_irad->GetYaxis()->SetTitle("#sigma(reco - gen) / gen");
  outFile->cd();
  g_res_vs_irad->Write();

  TDirectory* dirLayer = outFile->mkdir("resVsLayer");
  dirLayer->cd();
  TGraphErrors* g_res_vs_layer = makeResGraph(h_enDiff_vs_layer, 8, 26, 1);
  g_res_vs_layer->SetName("g_res_vs_layer");
  g_res_vs_layer->GetXaxis()->SetTitle("layer");
  g_res_vs_layer->GetYaxis()->SetTitle("#sigma(reco - gen) / gen");
  outFile->cd();
  g_res_vs_layer->Write();

  for(auto itr=hMap_enDiff_vs_layer.begin(); itr!=hMap_enDiff_vs_layer.end(); ++itr)
  {
    TDirectory* dir = outFile->mkdir(Form("resVsLayer_%f",itr->first));
    dir->cd();
    g_res_vs_layer = makeResGraph(itr->second, 8, 26, 1);
    g_res_vs_layer->SetName(Form("g_res_vs_layer_%f",itr->first));
    g_res_vs_layer->GetXaxis()->SetTitle("layer");
    g_res_vs_layer->GetYaxis()->SetTitle("#sigma(reco - gen) / gen");
    outFile->cd();
    g_res_vs_layer->Write();
  }


  //theApp -> Run();
  outFile->Write();
  outFile->Close();
  return 0;
}
