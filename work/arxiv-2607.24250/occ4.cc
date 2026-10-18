
#include "TApplication.h"

#include "iostream"
#include "iomanip"

#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"

#include "PMTStyle.h"
#include "PMType.h"

#include "Pedestal.h"
#include "SPEResponse.h"
#include "PMT.h"

#include "Occupancy.h"

using namespace std;

Double_t wbin;

Double_t  my_g( Double_t *x, Double_t *par)
{
  Double_t xx = x[0];

  Double_t Norm = par[0];
  Double_t Q0 = par[1];
  Double_t s0 = par[2];
        
  Double_t arg = 0.0; 
  if ( s0!=0 ) arg = ( xx - Q0 )/s0;    
  else cout << "Error: The code tries to divide by zero." << endl;
  
  Double_t result =  wbin*Norm/( sqrt( 2.0*TMath::Pi() )*s0 ) * TMath::Exp( -0.5*arg*arg );
  
  return result;

}

Int_t occ2( Float_t mu )
{
  time_t start;  
  time( &start );
  cout << "" << endl;
  cout << " The macro starts ( occ2.C ) ... " << endl;
  cout << "" << endl;
  gROOT->Reset();
  
  PMTStyle::SetDefaultStyle();

  
  //TCanvas *c1 = new TCanvas( "c1", "" );
  //c1->cd();
  //c1->SetLogy();

  Double_t Q0 = 0.02199;
  Double_t s0 = 0.00562;
  Pedestal ped( Q0, s0 );
  
  Double_t lambda = 7.79;
  Double_t theta = 5.1;
  
  Double_t alpha = 23;
  Double_t w = 0.41;

  Double_t alpha1 = 46;
  Double_t w1 = 0.20;
    
  Double_t p[6] = { lambda, theta, alpha, w, alpha1, w1 };
  SPEResponse gamma( PMType::GAMMA2EXP, p );
  
  Int_t nbins = 1000;
  Double_t xmin = -0.6395;
  Double_t xmax = 1.3605;

  //my_func->SetParNames( "Norm", "#lambda", "#theta" );

  TF1 *my_func  = new TF1( "my_func" , my_g, xmin, xmax, 3 );
  my_func->SetParNames( "Norm", "Q0", "s0" );
  my_func->SetLineColor( kRed );
  my_func->SetNpx( nbins );


  TFile *f1 = new TFile( Form( "out/test-occ4-%0.2f.root", mu ), "RECREATE" );
  f1->cd();

  TH1D *hG = new TH1D( "hG", "", 400, -50, 50 );
  TH1D *hmu = new TH1D( "hmu", "", 400, -50, 50 );

  hG->Sumw2();
  hmu->Sumw2();
  
  Int_t ntoys = 1000;
  for ( Int_t j=0; j<ntoys; j++ )
    {
  PMT specimen0( nbins, xmin, xmax, ped, gamma );
  
  specimen0.GenSpectrum( 1.0e+4, 0.0 );
  specimen0.GetSpectrum()->SetStats(0);
  //specimen0.DrawSpectrum();
    
  wbin = specimen0.GetSpectrum()->GetXaxis()->GetBinWidth(1);

  my_func->SetParameter( 0, specimen0.GetSpectrum()->GetEntries() ); 
  my_func->SetParameter( 1, specimen0.GetSpectrum()->GetMean() );
  my_func->SetParameter( 2, specimen0.GetSpectrum()->GetRMS() ); 

  specimen0.GetSpectrum()->Fit( "my_func", "", "", xmin, xmax );
  Double_t *par = my_func->GetParameters();
  
  //c1->Update();
  //c1->WaitPrimitive();

  
  PMT specimen( nbins, xmin, xmax, ped, gamma );
  Int_t ntot = 2.5e+5;

  specimen.GenSpectrum( ntot, mu );
  //specimen.GetSpectrum()->SetStats(0);
  //specimen.DrawSpectrum();
  
  //c1->Update();
  //c1->WaitPrimitive();

  Double_t Gtrue = ( w/alpha + w1/alpha1 + (1.0-w-w1)/lambda );
    
  /* ... */
  
  Occupancy occ( par[1], par[2] );
  Float_t Gfit = occ.FindG( specimen.GetSpectrum(), 0.2 );
  
  /* ... */

  Float_t DG = ( Gfit/Gtrue - 1.0 )*100.0;
  Float_t Dmu = ( occ.mu_bf/mu - 1.0 )*100.0;

  cout << j << endl;
  cout << mu << endl;
  
  cout << " True Gain : " << Gtrue << endl;
  cout << " BF Gain   : " << Gfit  << endl;
  cout << " Deviation : " << DG    << endl;
  cout << " Deviation : " << Dmu   << endl;
  cout << "" << endl;

  hG->Fill( DG );
  hmu->Fill( Dmu );
  //getchar();
    }

  cout << hG->GetMean() << " +/- " << hG->GetMeanError() << endl;
  cout << "" << endl;
  
  f1->Write();
  f1->Close();
    
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

int main( int argc, char *argv[] )  
{
  argc++;
  Float_t mu = atof( argv[1] );
      
  TApplication theApp( "App", 0, 0 );
  Int_t status = occ2( mu );
  return status;
  
}

