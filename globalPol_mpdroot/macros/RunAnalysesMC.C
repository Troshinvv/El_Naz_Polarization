void RunAnalysesMC(TString output){

   gROOT->LoadMacro("mpdloadlibs.C");
   gROOT->ProcessLine("mpdloadlibs()");

   MpdAnalysisManager man("ManagerAnal",-1) ;
   man.InputFileList("list.txt") ;
   man.ReadBranches("*") ; 
   
   MpdCentralityAll pCentr("pCentr","pCentr") ;
   man.AddTask(&pCentr) ;

   MpdEventPlaneAll pEP("pEP","pEP") ;
   man.AddTask(&pEP) ;
   
   MpdGlobalPolarizationMC pGlobalPol("pGlobalPolMC",output) ;
   man.AddTask(&pGlobalPol) ;

   man.Process() ;

}
