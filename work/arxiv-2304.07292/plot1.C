
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
  //c1->SetGridx();
  //c1->SetGridy();
  

  const Int_t N = 6;
  
  Double_t x[N] = { 0.5, 1, 2, 3, 4, 5 }; //, 6 };
  
  Double_t y1[N] = { -0.01061,
		     -0.01445,
		     -0.06637,
		     -0.09986,
		     -0.1687,
		     -0.3445
		     //-0.3311,
  }; 

  Double_t er_y1[N] = { 0.01599,
			0.01445,
			0.0149,
			0.01892,
			0.03805,
			0.06561
			//0.07151
  };
  
  TGraphErrors *gr1 = new TGraphErrors( N, x, y1, 0, er_y1 );
  gr1->GetXaxis()->SetLimits( 0.0, 6.0 );
  gr1->SetMaximum( +8.9 );
  gr1->SetMinimum( -4.0 );
  gr1->SetMarkerStyle( 20 );
  gr1->SetMarkerSize( 0.8 );
  gr1->SetTitle( "; #mu; #Delta Q_{s}  in %;" );
  gr1->Draw( "APEZC" );
  
  c1->Update();
  c1->WaitPrimitive();

  const Int_t M = 6;
  Double_t y2[M] = { 0.1967,
		     0.2974,
		     0.6579,
		     1.143,
		     2.004,
		     3.7
		     
  }; 

  Double_t er_y2[M] = { 0.01318,
			0.01047,
			0.01844,
			0.04858,
			0.04945,
			0.0908
  };

  TGraphErrors *gr2 = new TGraphErrors( M, x, y2, 0, er_y2 );
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
  
  //c1->Update();
  //c1->WaitPrimitive();
  
  
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
