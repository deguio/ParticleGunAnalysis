
void compareResolution()
{

  //=== list scenarios to compare and assign colors
  std::map<std::string, std::string> labelsMap;
  std::map<std::string, int> colorsMap;

  std::string pathName =  "out/20190405/";
  std::string label = "";

  label = "caliceDigi";
  labelsMap[label] = pathName+"PGunAnalysis_noiseScenario_0_algo_1_scaleByArea_False.root";
  colorsMap[label] = kGreen;

  label = "newDigi";
  labelsMap[label] = pathName+"PGunAnalysis_noiseScenario_0_algo_2_scaleByArea_False.root";
  colorsMap[label] = kBlue;

  label = "newDigi scaled (geometry)";
  labelsMap[label] = pathName+"PGunAnalysis_noiseScenario_0_algo_2_scaleByArea_True.root";
  colorsMap[label] = kRed;

  label = "newDigi scaled (dose)";
  labelsMap[label] = pathName+"PGunAnalysis_noiseScenario_3000_algo_2_scaleByArea_False.root";
  colorsMap[label] = kBlue+2;

  label = "newDigi scaled (geometry+dose)";
  labelsMap[label] = pathName+"PGunAnalysis_noiseScenario_3000_algo_2_scaleByArea_True.root";
  colorsMap[label] = kRed+2;




  //=== list plots to compare for each scenario and label axis
  //                     plotName,   xAxisTitle, yAxisTitle
  std::vector<std::tuple<std::string,std::string,std::string>> plotsVector;
  plotsVector.push_back(std::make_tuple("g_res_vs_irad","IRadius","#sigma((reco - gen) / gen)"));
  plotsVector.push_back(std::make_tuple("g_res_vs_layer","Layer","#sigma((reco - gen) / gen)"));

  plotsVector.push_back(std::make_tuple("g_res_vs_layer_0.000000","Layer (iRad<=10)","#sigma((reco - gen) / gen)"));
  plotsVector.push_back(std::make_tuple("g_res_vs_layer_10.000000","Layer (iRad>10 && iRad<=30)","#sigma((reco - gen) / gen)"));
  plotsVector.push_back(std::make_tuple("g_res_vs_layer_30.000000","Layer (iRad>30)","#sigma((reco - gen) / gen)"));

  for(auto plotsItr=plotsVector.begin(); plotsItr!=plotsVector.end(); ++plotsItr)
  {

    TCanvas* cResVsIRad = new TCanvas();
    cResVsIRad->SetGridx();
    cResVsIRad->SetGridy();
    TMultiGraph* mgResVsIRad = new TMultiGraph();
    TLegend* leg = new TLegend(0.2, 0.7, 0.7, 0.9);

    int count=0;
    for(auto labelsItr=labelsMap.begin(); labelsItr!=labelsMap.end(); ++labelsItr)
    {
      std::string key = labelsItr->first;
      TFile* file = TFile::Open((labelsItr->second).c_str());
      TGraphErrors* gResVsIRad = (TGraphErrors*)file->Get(std::get<0>(*plotsItr).c_str());
      gResVsIRad->SetTitle(key.c_str());
      gResVsIRad->SetMarkerColor(colorsMap[key]);
      gResVsIRad->SetLineColor(colorsMap[key]);
      mgResVsIRad->Add(gResVsIRad);
      leg->AddEntry(gResVsIRad, key.c_str());

      ++count;
    }

    mgResVsIRad->GetXaxis()->SetTitle(std::get<1>(*plotsItr).c_str());
    mgResVsIRad->GetYaxis()->SetTitle(std::get<2>(*plotsItr).c_str());
    mgResVsIRad->SetMinimum(0.1);
    mgResVsIRad->SetMaximum(0.3);
    mgResVsIRad->Draw("APL");
    leg->Draw("sames");
  }
  return;
}
