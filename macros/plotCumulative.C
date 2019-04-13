// .x macros/plotCumulative.C("HGCDigiIndex==2 && HGCDigiEta < 0")

void plotCumulative(std::string selection)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  tdrStyle->SetPadTopMargin(0.1);
  tdrStyle->SetPadTickX(0);

  //selection = "HGCDigiIndex==2 && HGCDigiEta < 0"; //consider only noise

  TLegend* leg = new TLegend(0.4, 0.75, 0.75, 0.85);


  TFile* fNew = TFile::Open("/Users/deguio/Desktop/HGCAL/Analysis/ntuples.nosync/20190413/noiseScenario_0_algo_2_scaleByArea_True_40.root");
  TTree* tNew = (TTree*)fNew->Get("hgcalTupleTree/tree");
  int nEvNew = tNew->GetEntries();

  TFile* fAged = TFile::Open("/Users/deguio/Desktop/HGCAL/Analysis/ntuples.nosync/20190413/noiseScenario_3000_algo_2_scaleByArea_True_117.root");
  TTree* tAged = (TTree*)fAged->Get("hgcalTupleTree/tree");
  int nEvAged = tAged->GetEntries();

  TH1F* hNew = new TH1F("hNew","hNew", 500, 0, 500);
  tNew->Draw("HGCDigiADC >> hNew", selection.c_str(), "goff");
  hNew->Scale(2./nEvNew); //consider 2 endcaps

  TH1F* hNewCumulative = (TH1F*)hNew->GetCumulative(kFALSE);
  hNewCumulative->GetXaxis()->SetTitle("Thr [ADC]");
  hNewCumulative->GetYaxis()->SetTitle("Channels above thr");
  hNewCumulative->SetLineWidth(3);
  hNewCumulative->SetLineColor(kBlue);
  leg->AddEntry(hNewCumulative, "new detector","L");


  TH1F* hAged = new TH1F("hAged","hAged", 500, 0, 500);
  tAged->Draw("HGCDigiADC >> hAged", selection.c_str(), "goff");
  hAged->Scale(2./nEvAged); //consider 2 endcaps

  TH1F* hAgedCumulative = (TH1F*)hAged->GetCumulative(kFALSE);
  hAgedCumulative->GetXaxis()->SetTitle("Thr [ADC]");
  hAgedCumulative->GetYaxis()->SetTitle("Channels above thr");
  hAgedCumulative->SetLineWidth(3);
  hAgedCumulative->SetLineColor(kRed+3);
  leg->AddEntry(hAgedCumulative, "aged detector [3000/fb, 16.7kh]","L");


  TCanvas* c1 = new TCanvas("c1","c1");
  c1->cd();
  c1->SetGridx();
  c1->SetGridy();
  c1->SetLogy();
  hAged->Draw("HISTO");
  hNew->Draw("HISTO,sames");

  TCanvas* c2 = new TCanvas("c2","c2");
  c2->cd();
  hAgedCumulative->Draw("HISTO");
  hNewCumulative->Draw("HISTO,sames");

  hAgedCumulative->GetXaxis()->SetRangeUser(0., 40.);
  hAgedCumulative->GetYaxis()->SetRangeUser(0.1,1e6);
  c2->Update();

  TGaxis *mipAxis = new TGaxis(gPad->GetUxmin(),gPad->GetUymax(),gPad->GetUxmax(),gPad->GetUymax(),0., 40./1024.*68.5, 510, "-L");
  mipAxis->SetLineColor(kRed);
  mipAxis->SetTextColor(kRed);
  mipAxis->SetLabelColor(kRed);
  mipAxis->SetTitle("Thr [MIP]");
  mipAxis->Draw("sames");
  leg->Draw("sames");

  c2->SetGridx();
  c2->SetGridy();
  c2->SetLogy();
  c2->Update();

  return;
}
