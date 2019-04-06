
void plotErrVsStat()
{
  int SoN = 0;
  int sig = 0;
  std::string label = "";

  std::vector<int> pointsVec = {500, 1000, 2500, 5000, 10000};
  std::string pathName =  "plotsBias_";


  std::map<std::string, std::pair<int, int>> labelsMap;
  std::map<std::string, int> colorsMap;

  SoN = 5; sig = 15;
  label = "SoN="+std::to_string(SoN)+", S="+std::to_string(sig);
  labelsMap[label] = std::make_pair(SoN, sig);
  colorsMap[label] = kRed;
  SoN = 5; sig = 30;
  label = "SoN="+std::to_string(SoN)+", S="+std::to_string(sig);
  labelsMap[label] = std::make_pair(SoN, sig);
  colorsMap[label] = kRed+2;
  SoN = 5; sig = 45;
  label = "SoN="+std::to_string(SoN)+", S="+std::to_string(sig);
  labelsMap[label] = std::make_pair(SoN, sig);
  colorsMap[label] = kRed+3;

//  SoN = 10; sig = 15;
//  label = "SoN="+std::to_string(SoN)+", S="+std::to_string(sig);
//  labelsMap[label] = std::make_pair(SoN, sig);
//  colorsMap[label] = kBlue;
//  SoN = 10; sig = 30;
//  label = "SoN="+std::to_string(SoN)+", S="+std::to_string(sig);
//  labelsMap[label] = std::make_pair(SoN, sig);
//  colorsMap[label] = kBlue+2;
//  SoN = 10; sig = 45;
//  label = "SoN="+std::to_string(SoN)+", S="+std::to_string(sig);
//  labelsMap[label] = std::make_pair(SoN, sig);
//  colorsMap[label] = kBlue+3;

  SoN = 20; sig = 15;
  label = "SoN="+std::to_string(SoN)+", S="+std::to_string(sig);
  labelsMap[label] = std::make_pair(SoN, sig);
  colorsMap[label] = kGreen+1;
  SoN = 20; sig = 30;
  label = "SoN="+std::to_string(SoN)+", S="+std::to_string(sig);
  labelsMap[label] = std::make_pair(SoN, sig);
  colorsMap[label] = kGreen+2;
  SoN = 20; sig = 45;
  label = "SoN="+std::to_string(SoN)+", S="+std::to_string(sig);
  labelsMap[label] = std::make_pair(SoN, sig);
  colorsMap[label] = kGreen+3;




  //=== list plots to compare for each scenario and label axis
  //                     plotName,   xAxisTitle, yAxisTitle
  std::vector<std::tuple<std::string,std::string,std::string>> plotsVector;
  plotsVector.push_back(std::make_tuple("biasErrorMap","n events","relative error on landau location par"));
  plotsVector.push_back(std::make_tuple("biasMap","n events","bias"));


  for(auto plotsItr=plotsVector.begin(); plotsItr!=plotsVector.end(); ++plotsItr)
  {
    TCanvas* cErrVsStat = new TCanvas();
    cErrVsStat->SetGridx();
    cErrVsStat->SetGridy();
    TMultiGraph* mgErrVsStats = new TMultiGraph();
    TLegend* leg = new TLegend(0.5, 0.7, 0.7, 0.9);

    for(auto labelsItr=labelsMap.begin(); labelsItr!=labelsMap.end(); ++labelsItr)
    {
      std::string key = labelsItr->first;

      TGraph* gErrVsStats = new TGraph();
      gErrVsStats->SetMarkerColor(colorsMap[key]);
      gErrVsStats->SetLineColor(colorsMap[key]);
      gErrVsStats->SetLineWidth(2);


      //loop over stats points
      int cont = 0;
      for(unsigned int pp=0; pp<pointsVec.size(); ++pp)
      {
        std::string statPoint = std::to_string(pointsVec[pp]);
        std::string fileName = pathName+statPoint+"/plots.root";
        TFile* file = TFile::Open(fileName.c_str());

        TH2F* hMap = (TH2F*)file->Get(std::get<0>(*plotsItr).c_str());

        int bin = hMap->FindBin(labelsItr->second.second, labelsItr->second.first);
        float val = hMap->GetBinContent(bin);
        gErrVsStats->SetPoint(cont, pointsVec[pp], val);

        ++cont;
      }

      mgErrVsStats->Add(gErrVsStats);
      std::string lab = "Sig = "+std::to_string(labelsItr->second.second)+", SoN = "+std::to_string(labelsItr->second.first);
      leg->AddEntry(gErrVsStats, lab.c_str());

    }


    mgErrVsStats->Draw("AL");
    mgErrVsStats->GetXaxis()->SetTitle(std::get<1>(*plotsItr).c_str());
    mgErrVsStats->GetYaxis()->SetTitle(std::get<2>(*plotsItr).c_str());
    leg->Draw("sames");

  }
  return;
}
