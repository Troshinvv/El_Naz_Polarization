void RunAnalysesRECO(TString output, TString analysis_choice, TString selection_choice){

   gROOT->LoadMacro("mpdloadlibs.C");
   gROOT->ProcessLine("mpdloadlibs()");

   MpdAnalysisManager man("ManagerAnal",1) ;
   man.InputFileList("list.txt") ;
   man.ReadBranches("*") ; 
   
   MpdCentralityAll pCentr("pCentr","pCentr") ;
   man.AddTask(&pCentr) ;

   MpdEventPlaneAll pEP("pEP","pEP") ;
   man.AddTask(&pEP) ;
   
   MpdGlobalPolarizationRECO pGlobalPol("pGlobalPolRECO",output,analysis_choice,selection_choice) ;
   man.AddTask(&pGlobalPol) ;

   man.Process() ;

}
