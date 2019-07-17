// root 'macros/plotCumulativeRDF.C("HGCDigiADC[HGCDigiIndex==2 && HGCDigiEta < 0]", 64)'

double xMin = 0.;
double xMax = 25.;
double yMin = 0.0001;
double yMax = 1.2;

std::map<int, std::vector<float>> hgcrocMap_;
std::map<int, std::vector<int>> hgcrocNcellsMap_;

std::string selection = "HGCDigiIndex==2 && HGCDigiEta < 0";


void createBinning()
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


void plotProbAboveNoiseRDF()
{
  ROOT::EnableImplicitMT();

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  tdrStyle->SetPadTopMargin(0.1);
  tdrStyle->SetPadTickX(0);
  tdrStyle->SetPadTickY(0);

  TFile* outFile = new TFile("out/plotProbAboveNoiseRDF.root","RECREATE");
  outFile->mkdir("plotter");
  outFile->cd("plotter");

  createBinning();


  ROOT::RDataFrame dfAged("hgcalTupleTree/tree","/Users/deguio/Desktop/HGCAL/Analysis/ntuples.nosync/20190511/noiseScenario_3000_algo_2_pileup_0_scaleByArea_True_4.root");
  auto nEvAged = dfAged.Filter("event>0").Count();

  int firstLayer = hgcrocMap_.begin()->first;
	int lastLayer = hgcrocMap_.rbegin()->first;
  std::string radius = "10*sqrt(HGCDigiPosx*HGCDigiPosx + HGCDigiPosy*HGCDigiPosy)"; //rad in mm
  std::string ped = "0.25*(HGCDigiADCSum - HGCDigiADC)"; //average over first and last 2 samples of the pulse

  //first loop: compute average pedestal from first 2 and last 2 samples
  std::map<int, ROOT::RDF::RResultPtr<TProfile>> pedestal_layerMap;
  for(int lay=firstLayer; lay<lastLayer+1; ++lay)
  {
    std::string hName = "pedestal_layer"+std::to_string(lay);
    std::string mySelection = "[" +selection + " && HGCDigiLayer==" + std::to_string(lay)+"]";

    std::string pedVal = ped + mySelection;
    std::string radVal = radius + mySelection;

    auto filtered = dfAged.Define("myRad", radVal.c_str())
                          .Define("myPed", pedVal.c_str());

    pedestal_layerMap[lay] = filtered.Profile1D({hName.c_str(),hName.c_str(), int(hgcrocMap_[lay].size()-1), hgcrocMap_[lay].data()}, "myRad", "myPed");
  }

  //here is where the fill of PED plots happens
  for(auto elem : pedestal_layerMap)
    elem.second->Write();


  //second loop: book per layer plots
  std::map<int, ROOT::RDF::RResultPtr<TProfile>> probNoiseAboveHalfMip_layerMap;
  for(int lay=firstLayer; lay<lastLayer+1; ++lay)
  {
    std::string hName = "probNoiseAboveHalfMip_layer"+std::to_string(lay);
    std::string mySelection = "[" +selection + " && HGCDigiLayer==" + std::to_string(lay)+"]";

    std::string radVal = radius + mySelection;
    std::string sigVal = "HGCDigiADC" + mySelection;

    auto filtered = dfAged.Define("myRad", radVal.c_str())
                          .Define("mySig", sigVal.c_str());
    probNoiseAboveHalfMip_layerMap[lay] = filtered.Profile1D({hName.c_str(),hName.c_str(), int(hgcrocMap_[lay].size()-1), hgcrocMap_[lay].data()}, "myRad", "mySig");
  }

  for(auto elem : probNoiseAboveHalfMip_layerMap)
  {
    elem.second->Write();

    elem.second->Add(pedestal_layerMap[elem.first].GetPtr(), -1.);

    for(int bin=1; bin<elem.second->GetXaxis()->GetNbins()+1; ++bin)
    {
      double binContent = elem.second->GetBinContent(bin);
      elem.second->SetBinContent(bin, binContent/36./hgcrocNcellsMap_[elem.first][bin-1]);
    }
    //elem.second->Write();
  }

  return;
}
