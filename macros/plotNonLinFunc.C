{

  //27370 pix
  TF1* f_corr_sim = new TF1("f_corr_sim", "pol2", 0, 50000);
  f_corr_sim->SetParameters(1.009850e+00, 8.075250e-06, 2.905460e-10);
  f_corr_sim->SetLineColor(kBlue);
  f_corr_sim->SetNpx(10000);

  TF1* f_corr_data = new TF1("f_corr_data", "pol2", 0, 50000);
  f_corr_data->SetParameters(1.000000, 2.71238e-05, 1.32877e-10);
  f_corr_data->SetLineColor(kRed);
  f_corr_data->SetNpx(10000);

  TCanvas* c27370 = new TCanvas("c27370","c27370");
  f_corr_sim->Draw();
  f_corr_data->Draw("sames");

}
