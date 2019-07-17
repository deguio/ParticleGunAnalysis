void plotDeDxPeak(std::string p)
{
  gStyle->SetOptFit();
  tdrStyle->SetPadRightMargin(0.26);


  TChain* tree = new TChain("hgcalTupleTree/tree");
  tree->Add(Form("/Users/deguio/Desktop/HGCAL/Analysis/ntuples.nosync/muGun_sim_%sGeV_NTU.root",p.c_str()));

  THStack *hs = new THStack("hs","");
  int min=60;
  int max=120;

  TH1F* dedx120 = new TH1F("dedx120","dedx120",400,0,800);
  dedx120->SetLineColor(kRed);
  dedx120->SetMarkerColor(kRed);
  dedx120->SetLineWidth(2);
  TF1* mylandau_120 = new TF1("mylandau_120","landau",min,max);
  mylandau_120->SetNpx(10000);
  mylandau_120->SetLineColor(kRed);

  TH1F* dedx200 = new TH1F("dedx200","dedx200",400,0,800);
  dedx200->SetLineColor(kBlue);
  dedx200->SetMarkerColor(kBlue);
  dedx200->SetLineWidth(2);
  TF1* mylandau_200 = new TF1("mylandau_200","landau",min,max);
  mylandau_200->SetNpx(10000);
  mylandau_200->SetLineColor(kBlue);

  TH1F* dedx300 = new TH1F("dedx300","dedx300",400,0,800);
  dedx300->SetLineWidth(2);
  TF1* mylandau_300 = new TF1("mylandau_300","landau",min,max);
  mylandau_300->SetNpx(10000);
  mylandau_300->SetLineColor(kBlack);

  TCanvas* c120 = new TCanvas();
  tree->Draw("HGCSimHitsIntEnergy*1e9/3.62/120 >> dedx120","HGCSimHitsIntIndex!=2 && HGCSimHitsIntType==0","");
  dedx120->Scale(1./dedx120->GetEntries());
  dedx120->Fit(mylandau_120,"R");
  c120->Update();
  auto stats120 = (TPaveStats*)dedx120->GetListOfFunctions()->FindObject("stats");
  hs->Add(dedx120);

  TCanvas* c200 = new TCanvas();
  tree->Draw("HGCSimHitsIntEnergy*1e9/3.62/200 >> dedx200","HGCSimHitsIntIndex!=2 && HGCSimHitsIntType==1","");
  dedx200->Scale(1./dedx120->GetEntries());
  dedx200->Fit(mylandau_200,"R");
  c200->Update();
  auto stats200 = (TPaveStats*)dedx200->GetListOfFunctions()->FindObject("stats");
  hs->Add(dedx200);

  TCanvas* c300 = new TCanvas();
  tree->Draw("HGCSimHitsIntEnergy*1e9/3.62/300 >> dedx300","HGCSimHitsIntIndex!=2 && HGCSimHitsIntType==2","");
  dedx300->Scale(1./dedx120->GetEntries());
  dedx300->Fit(mylandau_300,"R");
  c300->Update();
  auto stats300 = (TPaveStats*)dedx300->GetListOfFunctions()->FindObject("stats");
  hs->Add(dedx300);


  TCanvas* c1 = new TCanvas();
  hs->Draw("NOSTACK");
  hs->GetXaxis()->SetTitle("dE/dx [e^{-}/#mum]");
  hs->GetXaxis()->SetRangeUser(0,250);
  hs->SetMaximum(0.07);
  c1->Modified();
  c1->Update();

  stats120->SetTextColor(kRed);
  stats200->SetTextColor(kBlue);
  stats300->SetTextColor(kBlack);
  stats120->SetX1NDC(0.77); stats120->SetX2NDC(0.97); stats120->SetY1NDC(0.65); stats120->SetY2NDC(0.9);
  stats200->SetX1NDC(0.77); stats200->SetX2NDC(0.97); stats200->SetY1NDC(0.40); stats200->SetY2NDC(0.65);
  stats300->SetX1NDC(0.77); stats300->SetX2NDC(0.97); stats300->SetY1NDC(0.15); stats300->SetY2NDC(0.4);

  c1->Modified();
  c1->Update();
  c1->Print(Form("~/Desktop/dedx/%s.png",p.c_str()));
  c1->Print(Form("~/Desktop/dedx/%s.pdf",p.c_str()));

  return;
}
