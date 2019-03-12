#include "interface/SetTDRStyle.h"
#include "interface/Utils.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TApplication.h"



int main(int argc, char** argv)
{
  gSystem -> Load("CfgManager/lib/libCFGMan.so");
  setTDRStyle();

  CfgManager opts;
  opts.ParseConfigFile(argv[1]);

  std::string inputName = opts.GetOpt<std::string>("Input.inputFile");
  std::string outputPath = opts.GetOpt<std::string>("Input.outputFolder");


  std::string createFolder = "mkdir -p "+outputPath;
  system(createFolder.c_str());
  std::string outputFile = outputPath+"/PGunAnalysis.root";
  TFile* outFile = new TFile(outputFile.c_str(), "RECREATE");
  outFile -> cd();

  TChain* chain = LoadChain(inputName);
  MyTreeVars tt;
  InitTreeVars(chain,tt);

  int nEntries = chain -> GetEntries();
  int nProcessed = 0;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if(entry % 1000 == 0)
    std::cout << "reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    chain -> GetEntry(entry);




    ++nProcessed;
  }
  std::cout << std::endl;



  outFile->Write();
  outFile->Close();
  return 0;
}
