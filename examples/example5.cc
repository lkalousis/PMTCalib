
#include "TApplication.h"

#include "iostream"
#include "iomanip"

#include "TCanvas.h"

#include "PMTStyle.h"
#include "PMType.h"

#include "Pedestal.h"
#include "SPEResponse.h"
#include "PMT.h"

#include "DFTmethod.h"

#include "SPEFitter.h"

using namespace std;

Int_t example5()
{
  time_t start;  
  
  time( &start );
  
  cout << "" << endl;
  
  cout << " The macro starts ( example5.C ) ... " << endl;

  cout << "" << endl;

  gROOT->Reset();
  
  PMTStyle::SetDefaultStyle();

  
  TCanvas *c1 = new TCanvas( "c1", "" );
  c1->cd();
  c1->SetLogy();

  Double_t Q0 = 20.0;
  Double_t s0 = 2.0;
  Pedestal ped( Q0, s0 );
  
  Double_t Q = 40.0;
  Double_t s = 13.0;
  Double_t alpha = 1.0/8.0;
  Double_t w = 0.2;
  Double_t p[4] = { Q, s, alpha, w };
  SPEResponse gaus( PMType::GAUSS, p );

  Int_t nbins = 200;
  Double_t xmin =   0.0;
  Double_t xmax = 800.0;
    
  PMT specimen( nbins, xmin, xmax, ped, gaus );
  Double_t mu = 0.25;
  Int_t ntot = 2.0e+5;
  specimen.GenSpectrum( ntot, mu );
  specimen.GetSpectrum()->SetStats(0);
  specimen.DrawSpectrum();

  //c1->Update();
  //c1->WaitPrimitive();
  
  
  SPEFitter fit;

  //Double_t mu_test = fit.FindMu( specimen.GetSpectrum(), Q0, s0 );
  //Double_t g_test = fit.FindG( specimen.GetSpectrum(), Q0, mu_test );
  
  Double_t p_test[4] = { Q, s, alpha, w };
  //Double_t p_test[4] = { g_test, 0.3*g_test, 1.0/(0.2*g_test), 0.2 };
  
  SPEResponse gaus_test( PMType::GAUSS, p_test );
  NumIntegration num( 2.0*nbins, xmin, xmax, gaus_test );
  
  
  num.wbin = specimen.GetSpectrum()->GetBinWidth(1);
  
  num.Norm = ntot;
  
  num.Q0 = Q0;
  num.s0 = s0;
  
  num.mu = mu;
  
  /*
  fit.SetNummethod( num );
  fit.FitwNummethod( specimen.GetSpectrum() );
  
  num.Norm = fit.vals[0];
  
  num.Q0 = fit.vals[1];
  num.s0 = fit.vals[2];

  num.mu = fit.vals[3]; 
  
  Double_t p_fit[4] = { fit.vals[4], fit.vals[5], fit.vals[6], fit.vals[7] };
  num.spef.SetParams( p_fit );
  */
  TGraph *grBF = num.GetGraph();
  grBF->Draw( "SAME,L" );
  /*
  TGraph *grPE[25];
    
  for ( Int_t i=0; i<25; i++ )
    {
      grPE[i] = num.GetGraphN( i );
      grPE[i]->Draw( "SAME,L" );
      
    }
  */
  Double_t Gtrue = ( w*alpha+(1.0-w)*Q );
  Double_t Gfit = ( fit.vals[7]*fit.vals[6]+(1.0-fit.vals[7])*fit.vals[4] ); 
  
  cout << " True Gain : " << Gtrue << endl;
  cout << " BF Gain   : " << Gfit  << endl;
  cout << " Deviation : " << ( Gfit/Gtrue - 1.0 )*100.0 << endl;
  
  cout << "" << endl;
  cout << "" << endl;

  while ( 1!=0 )
      {
  c1->Update();
  c1->WaitPrimitive();
      }
  
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

int main() 
{
  TApplication theApp( "App", 0, 0 );
  
  Int_t status = example5();

  return status;
  
}
