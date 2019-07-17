int nBins = 45000;
int rangeMax = 1e5;


std::pair<std::vector<float>,std::vector<float>> fitDeDxPeak(std::string p)
{
  std::vector<float> peaks;
  std::vector<float> means;
  TChain* tree = new TChain("hgcalTupleTree/tree");
  tree->Add(Form("/Users/deguio/Desktop/HGCAL/Analysis/ntuples.nosync/CloseByPGun/muGun_NTU_p%s.root",p.c_str()));

  THStack *hs = new THStack("hs","");
  int min=60;
  int max=120;

  TH1F* dedx120 = new TH1F("dedx120","dedx120",nBins,0,rangeMax);
  dedx120->SetLineColor(kRed);
  dedx120->SetMarkerColor(kRed);
  dedx120->SetLineWidth(2);
  TF1* mylandau_120 = new TF1("mylandau_120","landau",min,max);
  mylandau_120->SetNpx(10000);
  mylandau_120->SetLineColor(kRed);

  TH1F* dedx200 = new TH1F("dedx200","dedx200",nBins,0,rangeMax);
  dedx200->SetLineColor(kBlue);
  dedx200->SetMarkerColor(kBlue);
  dedx200->SetLineWidth(2);
  TF1* mylandau_200 = new TF1("mylandau_200","landau",min,max);
  mylandau_200->SetNpx(10000);
  mylandau_200->SetLineColor(kBlue);

  TH1F* dedx300 = new TH1F("dedx300","dedx300",nBins,0,rangeMax);
  dedx300->SetLineWidth(2);
  TF1* mylandau_300 = new TF1("mylandau_300","landau",min,max);
  mylandau_300->SetNpx(10000);
  mylandau_300->SetLineColor(kBlack);

  TCanvas* c120 = new TCanvas();
  tree->Draw("HGCSimHitsIntEnergy*1e9/3.62/120 >> dedx120","HGCSimHitsIntLayer<11 && HGCSimHitsIntIndex!=2 && HGCSimHitsIntType==0","goff");
  dedx120->Scale(1./dedx120->GetEntries());
  dedx120->Fit(mylandau_120,"QR");
  c120->Update();
  auto stats120 = (TPaveStats*)dedx120->GetListOfFunctions()->FindObject("stats");
  hs->Add(dedx120);

  TCanvas* c200 = new TCanvas();
  tree->Draw("HGCSimHitsIntEnergy*1e9/3.62/200 >> dedx200","HGCSimHitsIntLayer<11 && HGCSimHitsIntIndex!=2 && HGCSimHitsIntType==1","goff");
  dedx200->Scale(1./dedx120->GetEntries());
  dedx200->Fit(mylandau_200,"QR");
  c200->Update();
  auto stats200 = (TPaveStats*)dedx200->GetListOfFunctions()->FindObject("stats");
  hs->Add(dedx200);

  TCanvas* c300 = new TCanvas();
  tree->Draw("HGCSimHitsIntEnergy*1e9/3.62/300 >> dedx300","HGCSimHitsIntLayer<11 && HGCSimHitsIntIndex!=2 && HGCSimHitsIntType==2","goff");
  dedx300->Scale(1./dedx120->GetEntries());
  dedx300->Fit(mylandau_300,"QR");
  c300->Update();
  auto stats300 = (TPaveStats*)dedx300->GetListOfFunctions()->FindObject("stats");
  hs->Add(dedx300);


  TCanvas* c1 = new TCanvas();
  hs->Draw("NOSTACK");
  hs->GetXaxis()->SetTitle("dE/dx [e^{-}/#mum]");
  hs->GetXaxis()->SetRangeUser(0,250);
  hs->SetMaximum(0.1);
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

  means.push_back(dedx120->GetMean());
  means.push_back(dedx200->GetMean());
  means.push_back(dedx300->GetMean());
  peaks.push_back(mylandau_120->GetParameter(1));
  peaks.push_back(mylandau_200->GetParameter(1));
  peaks.push_back(mylandau_300->GetParameter(1));

  delete dedx120;
  delete dedx200;
  delete dedx300;
  delete c120;
  delete c200;
  delete c300;
  return std::make_pair(means,peaks);
}



void plotBB()
{
  gStyle->SetOptFit();
  tdrStyle->SetPadRightMargin(0.26);

  float muMass = 0.1056;
  //std::vector<std::string> momentum = {"0.4","0.6"};
  std::vector<std::string> momentum = {"0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","4.0","6.0","8.0","10.0","20.0","40.0","60.0","80.0","100.0","150.0","200.0"};

  TLegend* legBB = new TLegend(0.82, 0.6, 0.95, 0.85);
  TLegend* legRatBB = new TLegend(0.82, 0.6, 0.95, 0.85);

  TGraph* gBB120 = new TGraph();
  gBB120->SetMarkerColor(kRed);
  gBB120->SetLineColor(kRed);
  gBB120->SetLineWidth(2);
  TGraph* gBB200 = new TGraph();
  gBB200->SetMarkerColor(kBlue);
  gBB200->SetLineColor(kBlue);
  gBB200->SetLineWidth(2);
  TGraph* gBB300 = new TGraph();
  gBB300->SetLineWidth(2);

  TGraph* gBBmean120 = new TGraph();
  gBBmean120->SetMarkerColor(kRed);
  gBBmean120->SetLineColor(kRed);
  gBBmean120->SetLineWidth(2);
  gBBmean120->SetLineStyle(2);
  TGraph* gBBmean200 = new TGraph();
  gBBmean200->SetMarkerColor(kBlue);
  gBBmean200->SetLineColor(kBlue);
  gBBmean200->SetLineStyle(2);
  gBBmean200->SetLineWidth(2);
  TGraph* gBBmean300 = new TGraph();
  gBBmean300->SetLineWidth(2);
  gBBmean300->SetLineStyle(2);

  TGraph* gRatBB120 = new TGraph();
  gRatBB120->SetMarkerColor(kRed);
  gRatBB120->SetLineColor(kRed);
  gRatBB120->SetLineWidth(2);
  TGraph* gRatBB200 = new TGraph();
  gRatBB200->SetMarkerColor(kBlue);
  gRatBB200->SetLineColor(kBlue);
  gRatBB200->SetLineWidth(2);
  TGraph* gRatBB300 = new TGraph();
  gRatBB300->SetLineWidth(2);


  int count=0;
  for(auto p=momentum.begin(); p!=momentum.end(); ++p)
  {
    auto dedxs = fitDeDxPeak(*p);
    auto means = dedxs.first;
    auto peaks = dedxs.second;
    gBB120->SetPoint(count, std::stof(*p)/muMass, peaks[0]);
    gBB200->SetPoint(count, std::stof(*p)/muMass, peaks[1]);
    gBB300->SetPoint(count, std::stof(*p)/muMass, peaks[2]);

    gBBmean120->SetPoint(count, std::stof(*p)/muMass, means[0]);
    gBBmean200->SetPoint(count, std::stof(*p)/muMass, means[1]);
    gBBmean300->SetPoint(count, std::stof(*p)/muMass, means[2]);

    gRatBB120->SetPoint(count, std::stof(*p)/muMass, peaks[0]/107.18);
    gRatBB200->SetPoint(count, std::stof(*p)/muMass, peaks[1]/107.18);
    gRatBB300->SetPoint(count, std::stof(*p)/muMass, peaks[2]/107.18);

    ++count;
  }


  auto mgBB = new TMultiGraph();
  mgBB->Add(gBB120);
  mgBB->Add(gBB200);
  mgBB->Add(gBB300);
  mgBB->Add(gBBmean120);
  mgBB->Add(gBBmean200);
  mgBB->Add(gBBmean300);

  legBB->AddEntry(gBB120,"MPV 120#mum","L");
  legBB->AddEntry(gBB200,"MPV 200#mum","L");
  legBB->AddEntry(gBB300,"MPV 300#mum","L");
  legBB->AddEntry(gBBmean120,"Mean 120#mum","L");
  legBB->AddEntry(gBBmean200,"Mean 200#mum","L");
  legBB->AddEntry(gBBmean300,"Mean 300#mum","L");

  TCanvas* cBB = new TCanvas();
  cBB->SetLogx();
  cBB->SetGridx();
  cBB->SetGridy();
  mgBB->Draw("AL");
  mgBB->GetXaxis()->SetTitle("#beta#gamma");
  mgBB->GetYaxis()->SetTitle("dE/dx [e^{-}/#mum]");
  mgBB->GetXaxis()->SetRangeUser(2,1000);
  mgBB->SetMinimum(50);
  mgBB->SetMaximum(200);

  legBB->Draw("same");
  cBB->Modified();
  cBB->Update();
  cBB->Print("~/Desktop/dedx/bb.png");
  cBB->Print("~/Desktop/dedx/bb.pdf");


  //=== Ratio ===
  auto mgRatBB = new TMultiGraph();
  mgRatBB->Add(gRatBB120);
  mgRatBB->Add(gRatBB200);
  mgRatBB->Add(gRatBB300);

  legRatBB->AddEntry(gRatBB120,"120#mum","L");
  legRatBB->AddEntry(gRatBB200,"200#mum","L");
  legRatBB->AddEntry(gRatBB300,"300#mum","L");

  TCanvas* cRatBB = new TCanvas();
  cRatBB->SetLogx();
  cRatBB->SetGridx();
  cRatBB->SetGridy();
  mgRatBB->Draw("AL");
  mgRatBB->GetXaxis()->SetTitle("#beta#gamma");
  mgRatBB->GetYaxis()->SetTitle("MPV(dE/dx)/Mean(dE/dx)_{min}");
  mgRatBB->GetXaxis()->SetRangeUser(2,1000);
  mgRatBB->SetMinimum(0.5);
  mgRatBB->SetMaximum(1);

  legRatBB->Draw("same");
  cRatBB->Modified();
  cRatBB->Update();
  cRatBB->Print("~/Desktop/dedx/ratio_bb.png");
  cRatBB->Print("~/Desktop/dedx/ratio_bb.pdf");

  return;
}
