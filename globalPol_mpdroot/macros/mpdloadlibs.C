R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "gconfig/basiclibs.C"

void mpdloadlibs()
{
  /** Load basic and ROOT libraries **/
  basiclibs();

  /** Load FairRoot libraries **/
  // gSystem->Load("libBase");
  // gSystem->Load("libFairTools");
  // gSystem->Load("libGeoBase");
  // gSystem->Load("libParBase");
  // gSystem->Load("libGen");
  // gSystem->Load("libTrkBase");
  // gSystem->Load("libGeane");

  /** Load MpdRoot libraries **/
  // Base
  // gSystem->Load("libMpdBase");
  // gSystem->Load("libMpdMCStack");
  // gSystem->Load("libMpdField");
  // gSystem->Load("libPassive");
  // gSystem->Load("libUniGenFormat");
  // gSystem->Load("libMpdGen");
  // Hadgen
  // gSystem->Load("libHADGEN");
  // gSystem->Load("libTHadgen");
  gSystem->Load("libMpdGeneralGenerator");

  // Detectors
  // gSystem->Load("libtpc");
  // gSystem->Load("libTof");
  // gSystem->Load("libEtof"); // <-------------segmentation fault!!!!!!
  // gSystem->Load("libEmc");
  // gSystem->Load("libZdc");
  // gSystem->Load("libSts");
  // gSystem->Load("libFfd");
  //gSystem->Load("libFsa");
  //gSystem->Load("libBbc");
  //gSystem->Load("libNDet");
  //gSystem->Load("libStt");
  //gSystem->Load("libSft");

  // Reconstruction
  gSystem->Load("libKalman");
  gSystem->Load("libLHETrack");  // <----------- undefined symbol: _ZTI14TpcSectorGeoAZ
  gSystem->Load("libMpdPid");
  gSystem->Load("libMpdDst");
  gSystem->Load("libMpdMiniEvent");

  // MPD Physics
  gSystem->Load("libMpdPhysics.so");    // common
  gSystem->Load("libMpdMcDst");
  //gSystem->Load("libMpdFemtoMaker.so");      // femtoscopy
  //gSystem->Load("libMpdFemtoMakerUser.so");  // femtoscopy
}

// check whether file exists
bool CheckFileExist(TString& fileName)
{
    gSystem->ExpandPathName(fileName);
    if (gSystem->AccessPathName(fileName.Data()) == true)
    {
        cout<<endl<<"no specified file: "<<fileName<<endl;
        return false;
    }

    return true;
}
