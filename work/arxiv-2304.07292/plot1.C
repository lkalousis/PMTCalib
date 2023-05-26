
R__LOAD_LIBRARY(libPMTCalib)

#include "PMTStyle.h"

Int_t plot1()
{
  time_t start;  
  
  time( &start );
  
  cout << "" << endl;
  
  cout << " The macro starts ( plotter.C ) ... " << endl;

  cout << "" << endl;

  gROOT->Reset();
  PMTStyle::SetDefaultStyle();

  
  TCanvas *c1 = new TCanvas( "c1", "" );
  c1->cd();
  c1->SetGridx();
  c1->SetGridy();
  

  const Int_t N = 8;
  
  Double_t x[N] = { 0.5, 1, 2, 3, 4, 5, 6, 7 };//, 8 };
  
  Double_t y1[N] = { 0.00559,
		     -0.03223,
		     -0.04450,
		     -0.06647,
		     -0.09943,
		     -0.30910,
		     -0.30730,
		     -0.78270 //,
		     //-1.09600
  }; 

  Double_t er_y1[N] = { 0.01597,
			0.01447,
			0.01642,
			0.02017,
			0.02714,
			0.06258,
			0.07876,
			0.13360 //,
			//0.1876
  };
  
  TGraphErrors *gr1 = new TGraphErrors( N, x, y1, 0, er_y1 );
  gr1->GetXaxis()->SetLimits( 0.0, 8.0 );
  gr1->SetMaximum( +5.0 );
  gr1->SetMinimum( -5.0 );
  gr1->SetMarkerStyle( 20 );
  gr1->SetMarkerSize( 0.8 );
  gr1->SetTitle( "; #mu; #Delta Q_{s}  in %;" );
  gr1->Draw( "APEZC" );
  
  c1->Update();
  c1->WaitPrimitive();

  Double_t y2[6] = { 0.3458,
		     0.5008,
		     0.3077,
		     -0.2732,
		     -0.7901,
		     -1.488
		      
  }; 

  Double_t er_y2[6] = { 0.01474,
			0.01183,
			0.01350,
			0.01849,
			0.08911,
			0.11900
  };

  TGraphErrors *gr2 = new TGraphErrors( 6, x, y2, 0, er_y2 );
  gr2->SetMarkerStyle( 20 );
  gr2->SetMarkerSize( 0.8 );
  gr2->SetMarkerColor( 30 );
  gr2->SetLineColor( 30 );
  gr2->Draw( "PEZC" );

  TLegend *leg_1 = new TLegend( 0.2, 0.6, 0.5, 0.85 );
  leg_1->AddEntry( gr1, "DFT approach", "p" );
  leg_1->AddEntry( gr2, "Numerical integration", "p" );
  leg_1->SetFillColor( 0 ); 
  leg_1->Draw();
  
  c1->Update();
  c1->WaitPrimitive();

  
  cout << " ... the macro ends ! " << endl;
	  
  cout << "" << endl;
  
  time_t end;      
  
  time( &end );
      
  Int_t dura = difftime( end, start );      
   
  Int_t min = dura / 60; Int_t sec = dura % 60;
    
  cout << " ---> "<< Form( "%02d:%02d", min, sec ) << endl;  
    
  cout << "" << endl;

  cout << "" << endl;
    
  return 0;
  
}
