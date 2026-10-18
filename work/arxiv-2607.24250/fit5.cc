
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
#include "SPEFitter.h"
#include "PMT.h"

#include "Occupancy.h"
#include "TLatex.h"

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
  
  Double_t alpha1 = 32;
  Double_t w1 = 0.06;
    
  Double_t p[6] = { lambda, theta, alpha, w, alpha1, w1 };
  SPEResponse gamma( PMType::GAMMA2EXP, p );

  Int_t nbins = 1000;
  Double_t xmin = -0.6395;
  Double_t xmax = 1.3605;

  Double_t Gtrue = ( w/alpha + w1/alpha1 + (1.0-w-w1)/lambda );
  
  //my_func->SetParNames( "Norm", "#lambda", "#theta" );

  //TF1 *my_func  = new TF1( "my_func" , my_g, xmin, xmax, 3 );
  //my_func->SetParNames( "Norm", "Q0", "s0" );
  //my_func->SetLineColor( kRed );
  //my_func->SetNpx( nbins );


  TFile *f1 = new TFile( Form( "out/test-fit5-%0.2f.root", mu ), "RECREATE" );
  f1->cd();

  TH1D *hG = new TH1D( "hG", "", 400, -20, 20 );
  TH1D *hmu = new TH1D( "hmu", "", 400, -20, 20 );
  TH1D *hchi2 = new TH1D( "hchi2", "", 100, 0.0, 10.0 );
  
  hG->Sumw2();
  hmu->Sumw2();
  hchi2->Sumw2();
  
  Int_t ntoys = 1000;
  for ( Int_t j=0; j<ntoys; j++ )
    {
      /* PMT specimen0( nbins, xmin, xmax, ped, gamma );
  
  specimen0.GenSpectrum( 1.0e+4, 0.0 );
  specimen0.GetSpectrum()->SetStats(0);
  //specimen0.DrawSpectrum();
    
  wbin = specimen0.GetSpectrum()->GetXaxis()->GetBinWidth(1);

  my_func->SetParameter( 0, specimen0.GetSpectrum()->GetEntries() ); 
  my_func->SetParameter( 1, specimen0.GetSpectrum()->GetMean() );
  my_func->SetParameter( 2, specimen0.GetSpectrum()->GetRMS() ); 
      */
      //specimen0.GetSpectrum()->Fit( "my_func", "", "", xmin, xmax );
      //Double_t *parg = my_func->GetParameters();
  
  //c1->Update();
  //c1->WaitPrimitive();

  
  PMT specimen( nbins, xmin, xmax, ped, gamma );
  Int_t ntot = 2.5e+5;

  specimen.GenSpectrum( ntot, mu );
  specimen.GetSpectrum()->SetStats(0);
  
  Int_t bin = specimen.GetSpectrum()->GetMaximumBin();
  Double_t max = specimen.GetSpectrum()->GetBinContent( bin );
  specimen.GetSpectrum()->SetMinimum(0.1);
  specimen.GetSpectrum()->SetMaximum( max*2.5 );
  //specimen.GetSpectrum()->SetName( Form("h%d", j) );
  //specimen.DrawSpectrum();    
  //c1->Update();
  //c1->WaitPrimitive();
 
  /* ... */
  
  //Occupancy occ( par[1], par[2] );
  //Float_t Gfit = occ.FindG( specimen.GetSpectrum(), 0.2 );

  SPEFitter fit;
  Double_t mu_test = fit.FindMu( specimen.GetSpectrum(), Q0, s0 );
  Double_t g_test = fit.FindG( specimen.GetSpectrum(), Q0, mu_test );
  cout << mu_test << ", " << g_test << endl;

  Double_t p_test[4] = { 1.0/g_test, 7.0, 1.0/(0.5*g_test), 0.4 };
  SPEResponse gamma_test( PMType::GAMMA, p_test );
  DFTmethod dft( nbins, xmin, xmax, gamma_test );
  
  dft.wbin = specimen.GetSpectrum()->GetBinWidth(1);
  //cout << dft.wbin << endl;

  dft.Norm = specimen.GetSpectrum()->Integral();
  
  dft.Q0 = Q0;
  dft.s0 = s0;
  
  dft.mu = mu_test;
  //dft.CalculateValues();
  
  
  fit.SetDFTmethod( dft );
  fit.FitwDFTmethod( specimen.GetSpectrum() );
  
  Double_t Gfit = ( fit.vals[7]/fit.vals[6] + (1.0-fit.vals[7])/fit.vals[4] );
  
  dft.Norm = fit.vals[0];
  
  dft.Q0 = fit.vals[1];
  dft.s0 = fit.vals[2];

  dft.mu = fit.vals[3]; 
  
  //Double_t p_fit[4] = { fit.vals[4], fit.vals[5], fit.vals[6], fit.vals[7] };
  //dft.spef.SetParams( p_fit );
  /*
  TGraph *grBF = dft.GetGraph();
  grBF->Draw( "SAME,L" );

  TGraph *grPE[25];
    
  for ( Int_t i=0; i<25; i++ )
    {
      grPE[i] = dft.GetGraphN( i );
      grPE[i]->Draw( "SAME,L" );
      
    }
  
  TLatex *t = new TLatex();
  t->DrawLatex( 0.80, 1.0*max, Form( "#font[22]{#scale[0.8]{#mu = %.3f #pm %.3f}}", fit.vals[3], fit.errs[3] ) );
  TLatex *t1 = new TLatex();
  t1->DrawLatex( 0.80, 0.5*max, Form( "#font[22]{#scale[0.8]{#lambda = %.2f #pm %.2f}}", fit.vals[4], fit.errs[4] ) );
  TLatex *t2 = new TLatex();
  t2->DrawLatex( 0.80, 0.25*max, Form( "#font[22]{#scale[0.8]{#theta = %.1f #pm %.1f}}", fit.vals[5], fit.errs[5] ) );
  TLatex *t3 = new TLatex();
  t3->DrawLatex( 0.80, 0.13*max, Form( "#font[22]{#scale[0.8]{#alpha = %.0f #pm %.0f}}", fit.vals[6], fit.errs[6] ) );
  TLatex *t4 = new TLatex();
  t4->DrawLatex( 0.80, 0.075*max, Form( "#font[22]{#scale[0.8]{w = %.2f #pm %.2f}}", fit.vals[7], fit.errs[7] ) );
  TLatex *t5 = new TLatex();
  t5->DrawLatex( 0.80, 0.04*max, Form( "#font[22]{#scale[0.8]{#chi^{2}/NDOF = %.2f}}", fit.chi2r ) );
  */
  /* ... */

  Float_t DG = ( Gfit/Gtrue - 1.0 )*100.0;
  Float_t Dmu = ( fit.vals[3]/mu - 1.0 )*100.0;

  cout << j << endl;
  cout << " True Gain    : " << Gtrue << endl;
  cout << " BF Gain      : " << Gfit  << endl;
  cout << " Deviation DG : " << DG    << endl;
  cout << " Deviation mu : " << Dmu   << endl;
  cout << "" << endl;

  if ( fit.chi2r<3.0 )
	{
	  hG->Fill( DG );
	  hmu->Fill( Dmu );
	  hchi2->Fill( fit.chi2r );

	}

  
  //while (1!=0){
  //c1->Update();
  //c1->WaitPrimitive();
  //}
  
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

