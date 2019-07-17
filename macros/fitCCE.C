void fitCCE(std::string volts, std::string meas)
{
  std::string path = "cceData/"+volts+"V/cce_vs_f_300_"+volts+"V_"+meas+".txt";
  TGraph* mygraph_300 = new TGraph(path.c_str());
  path = "cceData/"+volts+"V/cce_vs_f_200_"+volts+"V_"+meas+".txt";
  TGraph* mygraph_200 = new TGraph(path.c_str());
  path = "cceData/"+volts+"V/cce_vs_f_120_"+volts+"V_"+meas+".txt";
  TGraph* mygraph_120 = new TGraph(path.c_str());

  TF1* pol1_300 = new TF1("pol1_300","pol1",0,2e16);
  TF1* pol1_200 = new TF1("pol1_200","pol1",0,2e16);
  TF1* pol1_120 = new TF1("pol1_120","pol1",0,2e16);

//  TF1* broken = new TF1("broken","x<[0] ? 1+[1]*x : [2]*x+([1]-[2])*[0]+1",0,1e17);
//  TF1* broken = new TF1("broken","x<=[0] ? 1+[2]*x : x>[0] && x<=[1] ? [3]*x+([2]-[3])*[0]+1 : [4]*x+([3]-[4])*[1]+([2]-[3])*[0]+1",0,1e17);
  TF1* broken = new TF1("broken","x<=[0] ? 1+[1]*x : (1-[2]*log(x))+([1]*[0]+[2]*log([0]))",0,1e17);

//  broken->FixParameter(0, 2.1e14);
//  broken->FixParameter(1, 7.5e14);


  TCanvas* c300 = new TCanvas("c300","c300");
  //c300.SetLogx();
  broken->FixParameter(0, 6e14);
  mygraph_300->Draw("AP");
  mygraph_300->Fit(broken,"R");
  //mygraph_300->GetXaxis()->SetLimits(1e14, 1e16);
  mygraph_300->GetYaxis()->SetRangeUser(0,1.2);

  TCanvas* c200 = new TCanvas("c200","c200");
  //c200.SetLogx();
  broken->FixParameter(0, 1.5e15);
//  broken->FixParameter(1, 1e17);
  mygraph_200->Draw("AP");
  mygraph_200->Fit(broken,"R");
  //mygraph_200->GetXaxis()->SetLimits(1e14, 1e16);
  mygraph_200->GetYaxis()->SetRangeUser(0,1.2);

  TCanvas* c120 = new TCanvas("c120","c120");
  //c120.SetLogx();
  broken->FixParameter(0, 1.5e15);
//  broken->FixParameter(1, 6e15);
  mygraph_120->Draw("AP");
  mygraph_120->Fit(broken,"R");
  //mygraph_120->GetXaxis()->SetLimits(1e14, 1e16);
  mygraph_120->GetYaxis()->SetRangeUser(0,1.2);
//
//  pol1_300->SetRange(1e14,1e16);
//  TCanvas fits;
//  pol1_300->Draw();
//  pol1_200->Draw("sames");
//  pol1_120->Draw("sames");



}
