void RunAnalyses(TString output, TString analysis_choice, TString selection_choice){

   gROOT->LoadMacro("mpdloadlibs.C");
   gROOT->ProcessLine("mpdloadlibs()");

   MpdAnalysisManager man("ManagerAnal",10) ;
   man.InputFileList("list.txt") ;
   man.ReadBranches("*") ; 
   
   MpdCentralityAll pCentr("pCentr","pCentr") ;
   man.AddTask(&pCentr) ;

   MpdEventPlaneAll pEP("pEP","pEP") ;
   man.AddTask(&pEP) ;
   
   //analysis_choice is either "selection" or "analysis" -> need to add the exit when you put something else! or change to enums!
   MpdGlobalPolarizationRECO pGlobalPol("pGlobalPolRECO",output,analysis_choice,selection_choice) ;
   //MpdGlobalPolarization_RECO pGlobalPol("pGlobalPol_RECO",output,analysis_choice,selection_choice) ;
   man.AddTask(&pGlobalPol) ;

   man.Process() ;

}
