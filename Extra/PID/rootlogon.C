/// Rootlogon for style

{
	TStyle* mcStyle = new TStyle("mcStyle","Liza's Root Styles"); 
	mcStyle->SetCanvasBorderMode(0);
	mcStyle->SetCanvasColor(kWhite);
	mcStyle->SetCanvasDefH(600); //Height of canvas
	mcStyle->SetCanvasDefW(600); //Width of canvas
	mcStyle->SetCanvasDefX(0);   //POsition on screen
	mcStyle->SetCanvasDefY(0);
	
	// For the Pad:
	mcStyle->SetPadBorderMode(0);
	mcStyle->SetPadColor(kWhite);
	mcStyle->SetPadGridX(false);
	mcStyle->SetPadGridY(false);
	mcStyle->SetGridColor(0);
	mcStyle->SetGridStyle(3);
	mcStyle->SetGridWidth(1);
  
	// For the frame:
	mcStyle->SetFrameBorderMode(0);
	mcStyle->SetFrameBorderSize(1);
	mcStyle->SetFrameFillColor(0);
	mcStyle->SetFrameFillStyle(0);
	mcStyle->SetFrameLineColor(1);
	mcStyle->SetFrameLineStyle(1);
	mcStyle->SetFrameLineWidth(1);
  
	// For the histo:
	mcStyle->SetHistLineColor(1);
	mcStyle->SetHistLineStyle(0);
	mcStyle->SetHistLineWidth(2);
	mcStyle->SetEndErrorSize(2);
    
	mcStyle->SetMarkerStyle(20);
  
	//For the fit/function:
//	mcStyle->SetOptFit(100);
	mcStyle->SetFitFormat("5.4g");
	mcStyle->SetFuncColor(2);
	mcStyle->SetFuncStyle(1);
	mcStyle->SetFuncWidth(1);

	//For the date:
	mcStyle->SetOptDate(0);

	// For the statistics box:
	mcStyle->SetOptFile(0);
//	mcStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
//	mcStyle->SetOptStat("mr");
//	mcStyle->SetStatColor(kWhite);
	mcStyle->SetStatColor(kGray);
	mcStyle->SetStatFont(42);
	mcStyle->SetStatFontSize(0.05);
	mcStyle->SetStatTextColor(1);
	mcStyle->SetStatFormat("6.4g");
	mcStyle->SetStatBorderSize(1);
	mcStyle->SetLegendBorderSize(1);
	mcStyle->SetStatColor(10); //white color
//	mcStyle->SetStatH(0.2);
//	mcStyle->SetStatW(0.25);
	mcStyle->SetStatX(0.97);
	mcStyle->SetStatY(0.98);
	mcStyle->SetStatH(0.03);
	mcStyle->SetStatW(0.3);

	// Margins:
	mcStyle->SetPadTopMargin(0.05);
	mcStyle->SetPadBottomMargin(0.125);
//	mcStyle->SetPadLeftMargin(0.1); // was original
	mcStyle->SetPadLeftMargin(0.125); // works better with current Title size
//	mcStyle->SetPadRightMargin(0.05);
	mcStyle->SetPadRightMargin(0.12); // this is for TH2D plots

	// For the Global title:
	mcStyle->SetOptTitle(0);
	mcStyle->SetTitleFont(42);
	mcStyle->SetTitleColor(1);
	mcStyle->SetTitleTextColor(1);
	mcStyle->SetTitleFillColor(10);
	mcStyle->SetTitleFontSize(0.05);

	// For the axis titles:
	mcStyle->SetTitleColor(1, "XYZ");
	mcStyle->SetTitleFont(42, "XYZ");
	mcStyle->SetTitleSize(0.05, "XYZ");
	mcStyle->SetTitleXOffset(0.85);
	mcStyle->SetTitleYOffset(1.0);
//	mcStyle->SetTitleYOffset(0.7);
	mcStyle->SetTitleSize(0.06,"xyz"); 

	// For the axis labels:
	mcStyle->SetLabelColor(1, "XYZ");
	mcStyle->SetLabelFont(42, "XYZ");
	mcStyle->SetLabelOffset(0.007, "XYZ");
//	mcStyle->SetLabelSize(0.04, "XYZ");
	mcStyle->SetLabelSize(0.05,"xyz");
	

	// For the axis:
	mcStyle->SetAxisColor(1, "XYZ");
	mcStyle->SetStripDecimals(kTRUE);
	mcStyle->SetTickLength(0.03, "XYZ");
//	mcStyle->SetNdivisions(510, "XYZ");
	mcStyle->SetNdivisions(505, "XYZ"); // control how many labels are written
//	mcStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
//	mcStyle->SetPadTickY(1);

	// Change for log plots:
	mcStyle->SetOptLogx(0);
	mcStyle->SetOptLogy(0);
	mcStyle->SetOptLogz(0);

	// Postscript options:
	mcStyle->SetPaperSize(20.,20.);
	mcStyle->SetPalette(55);
	//mcStyle->SetPalette(1,0); 
	mcStyle->SetHatchesLineWidth(5);
	mcStyle->SetHatchesSpacing(0.05);
	
	gROOT->SetStyle("mcStyle"); 
	cout << "Styles are Set!" << endl;     
	return;  
}
