// root 'macros/plotCumulativeRDF.C("HGCDigiADC[HGCDigiIndex==2 && HGCDigiEta < 0]")'

double xMin = 0.;
double xMax = 25.;
double yMin = 0.0001;
double yMax = 1.2;

void plotCumulativeRDF(std::string selection, int nChannels)
{
  ROOT::EnableImplicitMT();

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  tdrStyle->SetPadTopMargin(0.1);
  tdrStyle->SetPadTickX(0);
  tdrStyle->SetPadTickY(0);

  TLegend* leg = new TLegend(0.4, 0.75, 0.75, 0.85);


  ROOT::RDataFrame dfNew("hgcalTupleTree/tree","/Users/deguio/Desktop/HGCAL/Analysis/ntuples.nosync/20190416/noiseScenario_0_algo_2_scaleByArea_True_1.root");
  auto nEvNew = dfNew.Filter("event>0").Count();
  auto dfNew_filtered = dfNew.Define("myVar", selection.c_str());
  auto hNew = dfNew_filtered.Histo1D({"hNew","hNew",500,0,500}, "myVar");
  hNew->Scale(1./(*nEvNew)/36./nChannels);

  ROOT::RDataFrame dfAged("hgcalTupleTree/tree","/Users/deguio/Desktop/HGCAL/Analysis/ntuples.nosync/20190416/noiseScenario_3000_algo_2_scaleByArea_True_14.root");
  auto nEvAged = dfAged.Filter("event>0").Count();
  auto dfAged_filtered = dfAged.Define("myVar", selection.c_str());
  auto hAged = dfAged_filtered.Histo1D({"hAged","hAged",500,0,500}, "myVar");
  hAged->Scale(1./(*nEvAged)/36./nChannels);



  TH1F* hNewCumulative = (TH1F*)hNew->GetCumulative(kFALSE);
  hNewCumulative->GetXaxis()->SetTitle("Thr [ADC]");
  hNewCumulative->GetYaxis()->SetTitle("Ch fraction");
  hNewCumulative->SetLineWidth(3);
  hNewCumulative->SetLineColor(kBlue);
  leg->AddEntry(hNewCumulative, "new detector","L");

  TH1F* hAgedCumulative = (TH1F*)hAged->GetCumulative(kFALSE);
  hAgedCumulative->GetXaxis()->SetTitle("Thr [ADC]");
  hAgedCumulative->GetYaxis()->SetTitle("Ch fraction");
  hAgedCumulative->SetLineWidth(3);
  hAgedCumulative->SetLineColor(kRed+3);
  leg->AddEntry(hAgedCumulative, "aged detector [3000/fb, 16.7kh]","L");


  TCanvas* c2 = new TCanvas("c2","c2");
  c2->cd();
  hAgedCumulative->Draw("HISTO");
  hNewCumulative->Draw("HISTO,sames");

  hAgedCumulative->GetXaxis()->SetRangeUser(xMin, xMax);
  hAgedCumulative->GetYaxis()->SetRangeUser(yMin, yMax);
  c2->Update();

  TGaxis *mipAxis = new TGaxis(gPad->GetUxmin(),gPad->GetUymax(),gPad->GetUxmax(),gPad->GetUymax(), xMin, xMax/1024.*68.5, 510, "-L");
  mipAxis->SetLineColor(kRed);
  mipAxis->SetTextColor(kRed);
  mipAxis->SetLabelColor(kRed);
  mipAxis->SetTitle("Thr [MIP]");
  mipAxis->Draw("sames");

  TGaxis *chAxis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(), yMin*nChannels, yMax*nChannels, 510, "+G");
  chAxis->SetLineColor(kRed);
  chAxis->SetTextColor(kRed);
  chAxis->SetLabelColor(kRed);
  chAxis->SetLabelOffset(0.05);
  chAxis->SetTitleOffset(1.3);
  chAxis->SetTitle("nCh / HGROC");
  chAxis->Draw("sames");

  leg->Draw("sames");

  c2->SetGridx();
  c2->SetGridy();
  c2->SetLogy();
  c2->Update();

  std::string outPic = "occupancy/"+selection+".pdf";
  c2->Print(outPic.c_str());

  return;
}
