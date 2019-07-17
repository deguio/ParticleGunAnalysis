double yMin = 1e-6;
double yMax = 1.;
double xMin = 0;
double xMax = 25;


void plotProbAboveNoise_analytical(std::string infile)
{
  std::string folderName = "noiseProb_noise3000_pu0_fake";
  std::string command = "mkdir "+folderName;
  gSystem->Exec(command.c_str());
  TFile* inFile = new TFile(infile.c_str(),"READ");

  ofstream occupancyNumbers;
  occupancyNumbers.open(Form("%s/occupancyNumbers.txt", folderName.c_str()));

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend* leg = new TLegend(0.82, 0.3, 0.95, 0.85);

  TMultiGraph *mg = new TMultiGraph();
  std::map<int, THStack*> hs_map;
  std::map<int, TLegend*> leg_map;

  occupancyNumbers << std::setw(5) << "layerID"
                   << std::setw(15) << "rocID"
                   << std::setw(15) << "roc-radius"
                   << std::setw(15) << "chFrac/ev"
                   << std::setw(15) << "chN/ev\n";
  for (int lay=9; lay<23; ++lay)
  {
    std::string hname = "plotter/probNoiseAboveHalfMip_layer"+std::to_string(lay);
    TH1D* prof = (TH1D*)inFile->Get(hname.c_str());


    TGraphAsymmErrors* graph = new TGraphAsymmErrors();
    hs_map[lay] = new THStack(Form("hs_layer%d",lay),"");
    leg_map[lay] = new TLegend(0.82, 0.3, 0.95, 0.85);

    for(int roc=1; roc<prof->GetXaxis()->GetNbins()+1; ++roc)
    {
//      occupancyNumbers << std::fixed << std::setprecision(5) << std::setw(5) << lay
//                                                             << std::setw(15) << roc
//                                                             << std::setw(15) << prof->GetXaxis()->GetBinCenter(roc)
//                                                             << std::setw(15) << prof->GetBinContent(roc);
//
//
//      //plain occupancy
//      hname = "plotter/occRoc/occ_layer"+std::to_string(lay)+"_roc"+std::to_string(roc);
//      TH1D* hOcc_roc = (TH1D*)inFile->Get(hname.c_str());
//      occupancyNumbers << std::fixed << std::setprecision(5) << std::setw(15) << hOcc_roc->GetMean() << "\n";


      //compute prob
      graph->SetPoint(roc-1, prof->GetXaxis()->GetBinCenter(roc), prof->GetBinContent(roc));

//      //compute quantiles
//      const Int_t nq = 7;
//      Double_t xq[] = {0.01, 0.05, 0.32, 0.5, 0.68, 0.95, 0.99};
//      Double_t yq[nq];
//
//      hname = "plotter/probRoc/prob_layer"+std::to_string(lay)+"_roc"+std::to_string(roc);
//      TH1D* hProb_roc = (TH1D*)inFile->Get(hname.c_str());
//
//      hProb_roc->GetQuantiles(nq,yq,xq);
//      graph->	SetPointEYhigh(roc-1, yq[5]);
//
//
//      //compute cumulative distribution
//      hname = "plotter/sigRoc/sig_layer"+std::to_string(lay)+"_roc"+std::to_string(roc);
//      TH1D* hSig_roc = (TH1D*)inFile->Get(hname.c_str());
//      TH1D* hSig_rocCumulative = (TH1D*)hSig_roc->GetCumulative(kFALSE);
//      hSig_rocCumulative->SetLineColor(kOrange-2+2*roc);
//      hSig_rocCumulative->SetMarkerColor(kOrange-2+2*roc);
//      hSig_rocCumulative->SetMarkerSize(0.5);
//      hSig_rocCumulative->SetLineWidth(3);
//
//      hs_map[lay]->Add(hSig_rocCumulative);
//      leg_map[lay]->AddEntry(hSig_rocCumulative, Form("ROC %d",roc));

    }

    graph->SetLineColor(kOrange-12+lay);
    graph->SetMarkerColor(kOrange-12+lay);
    graph->SetMarkerSize(0.5);
    graph->SetLineWidth(3);

    graph->SetFillStyle(3005);
    graph->SetFillColor(kOrange-12+lay);

    mg->Add(graph);
    std::string lname = "layer "+std::to_string(lay);
    leg->AddEntry(graph,lname.c_str());

  }

  TCanvas* c1 = new TCanvas();
  c1->SetLogy();

  mg->GetYaxis()->SetTitle("fraction of channels above ZS per HGCROC");
  mg->GetXaxis()->SetTitle("rad [mm]");
  mg->SetMinimum(1e-6);
  mg->SetMaximum(0.1);

  mg->Draw("APL3");
  leg->Draw("L,same");

  c1->Print(Form("%s/noiseProb.pdf", folderName.c_str()));
  c1->Print(Form("%s/noiseProb.png", folderName.c_str()));

//  for(auto lay : hs_map)
//  {
//    TCanvas* c = new TCanvas(Form("c_lay%d",lay.first), Form("c_lay%d",lay.first));
//    c->cd();
//    gPad->SetTopMargin(0.1);
//    gPad->SetTickx(0);
//
//    lay.second->Draw("nostack");
//    lay.second->GetYaxis()->SetTitle("Fraction of channels");
//    lay.second->GetXaxis()->SetTitle("Thr [ADC]");
//    lay.second->GetXaxis()->SetRangeUser(xMin, xMax);
//    lay.second->SetMinimum(yMin);
//    lay.second->SetMaximum(yMax);
//    lay.second->Draw("nostack");
//    c->Update();
//    c->Modified();
//
//    TGaxis *mipAxis = new TGaxis(gPad->GetUxmin(),gPad->GetUymax(),gPad->GetUxmax(),gPad->GetUymax(), xMin, xMax/1024.*68.5, 510, "-L");
//    mipAxis->SetLineColor(kRed);
//    mipAxis->SetTextColor(kRed);
//    mipAxis->SetLabelColor(kRed);
//    mipAxis->SetTitle("Thr [MIP]");
//    mipAxis->Draw("sames");
//    c->Update();
//    c->Modified();
//
//    TLatex latex;
//    latex.SetTextSize(0.025);
//    latex.DrawLatexNDC(0.7, 0.85, Form("Layer %d",lay.first));
//    c->Update();
//    c->Modified();
//
//    leg_map[lay.first]->Draw("L,same");
//    c->Update();
//    c->Modified();
//
//    c->SetGridx();
//    c->SetGridy();
//    c->SetLogy();
//    c->Update();
//    c->Modified();
//
//    c->Print(Form("%s/cumulative_%d.pdf", folderName.c_str(), lay.first));
//    c->Print(Form("%s/cumulative_%d.png", folderName.c_str(), lay.first));
//
//  }

  return;
}
