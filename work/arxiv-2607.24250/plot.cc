
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

#include "DFTmethod.h"

using namespace std;

Int_t plot( Float_t mu )
{
  time_t start;  
  time( &start );
  cout << "" << endl;
  cout << " The macro starts ( plot.C ) ... " << endl;
  cout << "" << endl;
  gROOT->Reset();
  
  PMTStyle::SetDefaultStyle();

  
  TCanvas *c1 = new TCanvas( "c1", "" );
  c1->cd();
  //c1->SetLogy();
    
  Double_t Q0 = 0.02199;
  Double_t s0 = 0.00562;
  Pedestal ped( Q0, s0 );
  
  Double_t lambda = 7.79;
  Double_t theta = 5.1;
  
  Double_t alpha = 23;
  Double_t w = 0.41;
  
  Double_t p[4] = { lambda, theta, alpha, w };
  SPEResponse gamma( PMType::GAMMA, p );

  Int_t nbins = 1000;
  Double_t xmin = -0.14;
  Double_t xmax = 4000.0*0.0005+xmin;
  
  PMT specimen( nbins, xmin, xmax, ped, gamma );
  Int_t ntot = 2.5e+5;
  
  specimen.GenSpectrum( ntot, mu );
  specimen.GetSpectrum()->SetStats(0);
  specimen.DrawSpectrum();
  
  //c1->Update();
  //c1->WaitPrimitive();
  
  SPEResponse gamma_test( PMType::GAMMA, p );
  DFTmethod dft( 2.0*nbins, xmin, xmax, gamma_test );
      
  dft.wbin = specimen.GetSpectrum()->GetBinWidth(1);
  dft.Norm = ntot;
  
  dft.Q0 = ped.Q0;
  dft.s0 = ped.s0;
      
  dft.mu = mu;

  TGraph *grBF = dft.GetGraph();
  grBF->Draw( "AL" );

  Int_t npeaks = 10;
  TGraph *grPE[npeaks];
      
  for ( Int_t i=0; i<npeaks; i++ )
    {
      grPE[i] = dft.GetGraphN( i );
      grPE[i]->Draw( "SAME,L" );
      
    }
  while (1!=0){
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

int main( int argc, char *argv[] )  
{
  argc++;
  Float_t mu = atof( argv[1] );
    
  TApplication theApp( "App", 0, 0 );
  Int_t status = plot( mu );
  return status;
  
}

