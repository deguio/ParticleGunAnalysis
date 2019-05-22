{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TProfile2D* sonmap = (TProfile2D*)_file0->Get("plotter/signalToNoiseDoseAreaMap");
  TProfile2D* noisemap = (TProfile2D*)_file0->Get("plotter/noiseByFluenceMap");


  for (int bin=1; bin<sonmap->GetNumberOfBins()+1; ++bin)
  {
    if(sonmap->GetBinContent(bin) < 5.)
      sonmap->SetBinContent(bin, 0.0001);
  }

  TCanvas c1;
  sonmap->GetXaxis()->SetTitle("z [cm]");
  sonmap->GetYaxis()->SetTitle("r [cm]");
  sonmap->GetZaxis()->SetTitle("S/N");
  sonmap->Draw("COLZ");

  TCanvas c2;
  noisemap->GetXaxis()->SetTitle("z [cm]");
  noisemap->GetYaxis()->SetTitle("r [cm]");
  noisemap->GetZaxis()->SetTitle("Noise [PE]");
  noisemap->Draw("COLZ");


}
