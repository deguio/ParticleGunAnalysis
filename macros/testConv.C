{
  gSystem -> Load("lib/libPGunAnalysis.dylib");
  #include "interface/Utils.h"

  TF1 sigTimeModel("sigTimeModel", expoConv, 0, 25, 2);
  sigTimeModel.SetParameters(0.9,2.1);
  sigTimeModel.Draw();
}
