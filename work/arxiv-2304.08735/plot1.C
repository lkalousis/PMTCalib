
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
  
  Double_t y1[N] = { 0.107,
		     0.039,
		     0.009,
		     -0.040,
		     -0.003,
		     0.019
  }; 

  Double_t er_y1[N] = { 0.015,
			0.011,
			0.013,
			0.015,
			0.017,
			0.022
			//0.07151
  };
  
  TGraphErrors *gr1 = new TGraphErrors( N, x, y1, 0, er_y1 );
  gr1->GetXaxis()->SetLimits( 0.0, 6.0 );
  gr1->SetMaximum( +5.0 );
  gr1->SetMinimum( -3.0 );
  gr1->SetMarkerStyle( 20 );
  gr1->SetMarkerSize( 0.8 );
  gr1->SetTitle( "; #mu; #Delta Q_{s}  in %;" );
  gr1->Draw( "APEZC" );
  
  //c1->Update();
  //c1->WaitPrimitive();
  
  Double_t y2[N] = { 0.122,
		     0.020,
		     -0.030,
		     -0.121,
		     -0.163,
		     -0.120
		     
  }; 

  Double_t er_y2[N] = { 0.016,
			0.017,
			0.019,
			0.021,
			0.030,
			0.031
  };

  TGraphErrors *gr2 = new TGraphErrors( N, x, y2, 0, er_y2 );
  gr2->SetMarkerStyle( 20 );
  gr2->SetMarkerSize( 0.8 );
  gr2->SetMarkerColor( 30 );
  gr2->SetLineColor( 30 );
  gr2->Draw( "PEZC" );
  
  Double_t y3[N] = { 0.135,
		     -0.034,
		     -0.217,
		     -0.464,
		     -0.681,
		     -0.871
		     
  }; 

  Double_t er_y3[N] = { 0.015,
			0.013,
			0.015,
			0.019,
			0.029,
			0.041
  };

  TGraphErrors *gr3 = new TGraphErrors( N, x, y3, 0, er_y3 );
  gr3->SetMarkerStyle( 20 );
  gr3->SetMarkerSize( 0.8 );
  gr3->SetMarkerColor( 50 );
  gr3->SetLineColor( 50 );
  gr3->Draw( "PEZC" );
  
  
  TLegend *leg_1 = new TLegend( 0.2, 0.6, 0.5, 0.85 );
  leg_1->AddEntry( gr1, "#sigma/Q = 25 %", "p" );
  leg_1->AddEntry( gr2, "#sigma/Q = 35 %", "p" );
  leg_1->AddEntry( gr3, "#sigma/Q = 45 %", "p" );
  
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
