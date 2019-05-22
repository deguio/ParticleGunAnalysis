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

std::map<int, std::vector<float>> hgcrocMap_;
std::map<int, std::vector<int>> hgcrocNcellsMap_;
void createRocBinning();


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
  std::string noiseName = opts.GetOpt<std::string>("Input.noiseFile");
  std::string outputPath = opts.GetOpt<std::string>("Input.outputFolder");
  std::string outputName = opts.GetOpt<std::string>("Input.outputName");

  int detIndex = opts.GetOpt<int>("Params.detectorIndex");

  float gevPerMip = opts.GetOpt<float>("Params.gevPerMip");
  float mipRange = opts.GetOpt<float>("Params.mipRange");
  int nAdcBits = opts.GetOpt<int>("Params.nAdcBits");
  float zsThr = opts.GetOpt<float>("Params.zsThr");
  int nAdcChannels = pow(2, nAdcBits);

  float halfMip = nAdcChannels/mipRange/2.;

  std::string createFolder = "mkdir -p "+outputPath;
  system(createFolder.c_str());
  std::string outputFile = outputPath+"/"+outputName;
  TFile* outFile = new TFile(outputFile.c_str(), "RECREATE");
  outFile -> cd();

  TChain* chainNoise = LoadChain(noiseName);
  MyTreeVars tt;
  InitTreeVars(chainNoise,tt);

  //book histograms
  createRocBinning();
  int firstLayer = hgcrocMap_.begin()->first;
  int lastLayer = hgcrocMap_.rbegin()->first;

  //instantiate and initialize histograms
  std::map<int, TProfile*> ped_layerMap;
  std::map<int, TH1D*> prob_layerMap;
  std::map<int, TH1D*> occ_layerMap;

  std::map<int, std::map<int, TH1D*>> sig_rocMap;
  std::map<int, std::map<int, TH1D*>> occ_rocMap;
  std::map<int, std::map<int, TH1D*>> prob_rocMap;
  std::map<int, std::map<int, TProfile*>> ps_rocMap;

  std::map<int, std::map<int, double>> ped_channelMap;
  std::map<int, std::map<int, double>> count_rocMap;

  for(int lay=firstLayer; lay<lastLayer+1; ++lay)
  {
    ped_layerMap[lay]  = new TProfile(Form("pedestal_layer%d",lay),         "", int(hgcrocMap_[lay].size()-1), hgcrocMap_[lay].data());
    prob_layerMap[lay] = new TH1D(Form("probNoiseAboveHalfMip_layer%d",lay),"", int(hgcrocMap_[lay].size()-1), hgcrocMap_[lay].data());
    occ_layerMap[lay]  = new TH1D(Form("occ_layer%d",lay),                  "", int(hgcrocMap_[lay].size()-1), hgcrocMap_[lay].data());

    for(unsigned int roc=0; roc<hgcrocNcellsMap_[lay].size(); ++roc)
    {
      occ_rocMap[lay][roc+1] = new TH1D(Form("occ_layer%d_roc%d", lay, roc+1), "", 250, 0, 10);     //ADC
      prob_rocMap[lay][roc+1] = new TH1D(Form("prob_layer%d_roc%d", lay, roc+1), "", 100000, 0, 1);     //ADC
      sig_rocMap[lay][roc+1] = new TH1D(Form("sig_layer%d_roc%d", lay, roc+1), "", 50000, 0, 500);     //ADC
      ps_rocMap[lay][roc+1] = new TProfile(Form("ps_layer%d_roc%d", lay, roc+1), "", 5, 0, 5);     //ADC
    }
  }


  //FIRST LOOP ===================================================================================
  if(nEntries == -1)
    nEntries = chainNoise -> GetEntries();

  std::cout << "Compute PED" << std::endl;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if(entry % 10 == 0)
    std::cout << "reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    chainNoise -> GetEntry(entry);

    //loop over digis
    for(unsigned int digi=0; digi<tt.HGCDigiIEta->size(); ++digi)
    {
      int digiIndex = tt.HGCDigiIndex->at(digi);
      int digiIEta = tt.HGCDigiIEta->at(digi);

      //consider only the relevant subdet
      if(digiIndex != detIndex || digiIEta > 0.)
        continue;

      int digiLayer = tt.HGCDigiLayer->at(digi);
      float digiRad = 10 * sqrt(std::pow(tt.HGCDigiPosx->at(digi),2) + std::pow(tt.HGCDigiPosy->at(digi),2));
      float ped = (tt.HGCDigiSamples->at(digi)[0] + tt.HGCDigiSamples->at(digi)[1])/2;

      ped_layerMap[digiLayer]->Fill(digiRad, ped);              //per roc ped
      ped_channelMap[digiLayer][digiIEta] += ped/nEntries/288.; //per channel ped
    }
  }

  //SECOND LOOP ===================================================================================
  TChain* chain = LoadChain(inputName);
  InitTreeVars(chain,tt);

  if(nEntries == -1)
    nEntries = chain -> GetEntries();

  std::cout << std::endl;
  std::cout << "Run ANALYSIS" << std::endl;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if(entry % 10 == 0)
    std::cout << "reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    chain -> GetEntry(entry);

    //loop over digis
    for(unsigned int digi=0; digi<tt.HGCDigiIEta->size(); ++digi)
    {
      int digiIndex = tt.HGCDigiIndex->at(digi);
      int digiIEta = tt.HGCDigiIEta->at(digi);

      //consider only the relevant subdet
      if(digiIndex != detIndex || digiIEta > 0.)
        continue;

      int digiLayer = tt.HGCDigiLayer->at(digi);
      float digiRad = 10 * sqrt(std::pow(tt.HGCDigiPosx->at(digi),2) + std::pow(tt.HGCDigiPosy->at(digi),2));
      float sig = tt.HGCDigiSamples->at(digi)[2];

      int radBin = prob_layerMap[digiLayer]->GetXaxis()->FindBin(digiRad);
      //float ped = ped_layerMap[digiLayer]->GetBinContent(radBin);
      double ped = ped_channelMap[digiLayer][digiIEta];


      sig_rocMap[digiLayer][radBin]->Fill(sig-ped, 1./36./hgcrocNcellsMap_[digiLayer][radBin-1]/nEntries);

      if ((sig - ped) > halfMip)
      {
        occ_layerMap[digiLayer]->Fill(digiRad, 1./36./nEntries);
        prob_layerMap[digiLayer]->Fill(digiRad, 1./36./hgcrocNcellsMap_[digiLayer][radBin-1]/nEntries);

        count_rocMap[digiLayer][radBin] += 1./36;

        for(int ts=0; ts<5; ++ts)
          ps_rocMap[digiLayer][radBin]->Fill(ts, tt.HGCDigiSamples->at(digi)[ts]);

      }
      else
      {
        count_rocMap[digiLayer][radBin] += 0;
      }


    }//end loop digis

    //fill occupancy plots
    for(auto lay : count_rocMap)
      for(auto roc : lay.second)
      {
        occ_rocMap[lay.first][roc.first]->Fill(roc.second);
        prob_rocMap[lay.first][roc.first]->Fill(roc.second/hgcrocNcellsMap_[lay.first][roc.first-1]);

      }
    count_rocMap.clear();

  }

  std::cout << std::endl;


  outFile->mkdir("plotter/pedLayer");
  outFile -> cd("plotter/pedLayer");
  for(auto elem : ped_layerMap)
    elem.second->Write();

  outFile->mkdir("plotter/probLayer");
  outFile -> cd("plotter/probLayer");
  for(auto elem : prob_layerMap)
    elem.second->Write();

  outFile->mkdir("plotter/occLayer");
  outFile -> cd("plotter/occLayer");
  for(auto elem : occ_layerMap)
    elem.second->Write();

  outFile->mkdir("plotter/sigRoc");
  outFile -> cd("plotter/sigRoc");
  for(auto lay : sig_rocMap)
    for(auto roc : lay.second)
      roc.second->Write();

  outFile->mkdir("plotter/psRoc");
  outFile -> cd("plotter/psRoc");
  for(auto lay : ps_rocMap)
    for(auto roc : lay.second)
      roc.second->Write();

  outFile->mkdir("plotter/occRoc");
  outFile -> cd("plotter/occRoc");
  for(auto lay : occ_rocMap)
    for(auto roc : lay.second)
      roc.second->Write();

  outFile->mkdir("plotter/probRoc");
  outFile -> cd("plotter/probRoc");
  for(auto lay : prob_rocMap)
    for(auto roc : lay.second)
      roc.second->Write();



  outFile->cd();
  outFile->Close();
}



//=================================================================
void createRocBinning()
{
  std::vector<float> arr9 = {1537.0, 1790.7, 1997.1};
  hgcrocMap_[9] = arr9;
  std::vector<float> arr10 = {1537.0, 1790.7, 2086.2};
  hgcrocMap_[10] = arr10;
  std::vector<float> arr11 = {1537.0, 1790.7, 2132.2};
  hgcrocMap_[11] = arr11;
  std::vector<float> arr12 = {1537.0, 1790.7, 2179.2};
  hgcrocMap_[12] = arr12;
  std::vector<float> arr13 = {1378.2, 1503.9, 1790.7, 2132.2, 2326.6};
  hgcrocMap_[13] = arr13;
  std::vector<float> arr14 = {1378.2, 1503.9, 1790.7, 2132.2, 2430.4};
  hgcrocMap_[14] = arr14;
  std::vector<float> arr15 = {1183.0, 1503.9, 1790.7, 2132.2, 2538.8};
  hgcrocMap_[15] = arr15;
  std::vector<float> arr16 = {1183.0, 1503.9, 1790.7, 2132.2, 2594.8};
  hgcrocMap_[16] = arr16;
  std::vector<float> arr17 = {1183.0, 1503.9, 1790.7, 2132.2, 2594.8};
  hgcrocMap_[17] = arr17;
  std::vector<float> arr18 = {1183.0, 1503.9, 1790.7, 2132.2, 2594.8};
  hgcrocMap_[18] = arr18;
  std::vector<float> arr19 = {1037.8, 1157.5, 1503.9, 1790.7, 2132.2, 2594.8};
  hgcrocMap_[19] = arr19;
  std::vector<float> arr20 = {1037.8, 1157.5, 1503.9, 1790.7, 2132.2, 2594.8};
  hgcrocMap_[20] = arr20;
  std::vector<float> arr21 = {1037.8, 1157.5, 1503.9, 1790.7, 2132.2, 2594.8};
  hgcrocMap_[21] = arr21;
  std::vector<float> arr22 = {1037.8, 1157.5, 1503.9, 1790.7, 2132.2, 2484.0};
  hgcrocMap_[22] = arr22;
  std::vector<int> ncells9 = {64, 32};
  hgcrocNcellsMap_[9] = ncells9;
  std::vector<int> ncells10 = {64, 48};
  hgcrocNcellsMap_[10] = ncells10;
  std::vector<int> ncells11 = {64, 56};
  hgcrocNcellsMap_[11] = ncells11;
  std::vector<int> ncells12 = {64, 64};
  hgcrocNcellsMap_[12] = ncells12;
  std::vector<int> ncells13 = {40, 64, 64, 24};
  hgcrocNcellsMap_[13] = ncells13;
  std::vector<int> ncells14 = {40, 64, 64, 40};
  hgcrocNcellsMap_[14] = ncells14;
  std::vector<int> ncells15 = {88, 64, 64, 56};
  hgcrocNcellsMap_[15] = ncells15;
  std::vector<int> ncells16 = {88, 64, 64, 64};
  hgcrocNcellsMap_[16] = ncells16;
  std::vector<int> ncells17 = {88, 64, 64, 64};
  hgcrocNcellsMap_[17] = ncells17;
  std::vector<int> ncells18 = {88, 64, 64, 64};
  hgcrocNcellsMap_[18] = ncells18;
  std::vector<int> ncells19 = {40, 96, 64, 64, 64};
  hgcrocNcellsMap_[19] = ncells19;
  std::vector<int> ncells20 = {40, 96, 64, 64, 64};
  hgcrocNcellsMap_[20] = ncells20;
  std::vector<int> ncells21 = {40, 96, 64, 64, 64};
  hgcrocNcellsMap_[21] = ncells21;
  std::vector<int> ncells22 = {40, 96, 64, 64, 48};
  hgcrocNcellsMap_[22] = ncells22;
}
