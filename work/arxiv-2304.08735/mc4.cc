
#include "TApplication.h"

#include "iostream"
#include "iomanip"

#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"

#include "PMTStyle.h"
#include "PMType.h"

#include "Pedestal.h"
#include "SPEResponse.h"
#include "PMT.h"

#include "NumIntegration.h"

#include "SPEFitter.h"

using namespace std;

Int_t mc( Float_t mu )
{
  time_t start;  
  
  time( &start );
  
  cout << "" << endl;
  
  cout << " The macro starts ( mc.C ) ... " << endl;

  cout << "" << endl;

  gROOT->Reset();
  
  PMTStyle::SetDefaultStyle();

  
  TCanvas *c1 = new TCanvas( "c1", "" );
  c1->cd();
  c1->SetLogy();

  Double_t Q0 = 0.1500;
  Double_t s0 = 0.0028;
  Pedestal ped( Q0, s0 );
  
  Double_t Q = 0.030; Double_t per = 0.30;
  Double_t s = per*Q;
  Double_t alpha = 82.0;
  Double_t w = 0.2;
  Double_t p[4] = { Q, s, alpha, w };
  SPEResponse gaus( PMType::GAUSS, p );

  Int_t nbins = 250;
  Double_t xmin = 0.0;
  Double_t xmax = 1.0;
  
  //Double_t mu = 5.0;
  Int_t ntot = 2.5e+6;

  TFile *f = new TFile( Form( "./out/mc4_%.2f_%.2f.root", mu, per ), "RECREATE" );
  f->cd();
  
  TH1F *h_g = new TH1F( "h_g" , "; Gain in nVs; Entries", 100, 0.015, 0.035 );
  TH1F *h_dg = new TH1F( "h_dg" , "; G_{bf}/G_{true}-1 in %; Entries", 40, -20.0, 20.0 );
  TH1F *h_chi2 = new TH1F( "h_chi2" , "; #chi^{2}/NDOF; Entries",  150, 0.0, 15.0 );

  PMT specimen( nbins, xmin, xmax, ped, gaus );
  
  Int_t ntoys = 100;
  Int_t i=0;

  while ( i<ntoys )
  //for ( Int_t i=0; i<ntoys; i++ )
    {
      cout << " i   : " << i << endl;
      cout << " mu  : " << mu << endl;
      cout << " per : " << per << endl;
      cout << "" << endl; 
      
  specimen.GenSpectrum( ntot, mu );
  specimen.GetSpectrum()->SetName( Form( "h_%d", i ) );
  specimen.GetSpectrum()->SetStats(0);
  
  specimen.GetSpectrum()->GetXaxis()->SetRangeUser( 0.1, 0.7 );
  //specimen.GetSpectrum()->SetMaximum( 2.5e+5 );
  specimen.GetSpectrum()->SetMaximum( 2.5*specimen.GetSpectrum()->GetBinContent( specimen.GetSpectrum()->GetMaximumBin() ) );
  
  specimen.DrawSpectrum();
  
  
  SPEFitter fit;

  PMTModel mod( 2*nbins, xmin, xmax, PMType::ANATRUNCG );
  mod.wbin = specimen.GetSpectrum()->GetBinWidth(1);
  
  Double_t mu_test = fit.FindMu( specimen.GetSpectrum(), Q0, s0 );
  Double_t g_test = fit.FindG( specimen.GetSpectrum(), Q0, mu_test );
  
  //Double_t p_test[8] = { 1.0*ntot*1.0, Q0, s0, mu, Q, s, alpha, w };
  Double_t p_test[8] = { 1.0*ntot*1.0, Q0, s0, mu_test, g_test, 0.3*g_test, 1.0/(0.5*g_test), 0.2 };
  mod.SetParams( p_test );
  
  fit.SetPMTModel( mod );
  fit.FitwPMTModel( specimen.GetSpectrum() );
  
  Double_t p_bf[8] = { fit.vals[0], fit.vals[1], fit.vals[2], fit.vals[3], fit.vals[4], fit.vals[5], fit.vals[6], fit.vals[7] };
  mod.SetParams( p_bf );
  
  TGraph *grBF = mod.GetGraph();
  grBF->Draw( "SAME,L" );

  Int_t npeaks = 35;
  TGraph *grPE[35];
  
  for ( Int_t i=0; i<35; i++ )
    { //Int_t i = 10;
      grPE[i] = mod.GetGraphN( i );
      grPE[i]->Draw( "SAME,L" );
      
    } 

  TLatex *t = new TLatex();
  t->DrawLatex( 0.50, 1.0e+4, Form( "#font[22]{#scale[0.8]{#mu = %.3f #pm %.3f}}", fit.vals[3], fit.errs[3] ) );

  TLatex *t1 = new TLatex();
  t1->DrawLatex( 0.50, 0.5e+4, Form( "#font[22]{#scale[0.8]{#chi^{2}/NDOF = %.2f}}", fit.chi2r ) );

  
      
  Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );
  Double_t k = s/gn/sqrt( 2.0*TMath::Pi() )*TMath::Exp( -pow( Q, 2.0 )/( 2.0*pow( s, 2.0 ) ) );
  Double_t Qg = Q + k;
  Double_t G = ( w/alpha + (1.0-w)*Qg );
  
  
  Double_t Q_bf = fit.vals[4];
  Double_t s_bf = fit.vals[5];
  Double_t alpha_bf = fit.vals[6];
  Double_t w_bf = fit.vals[7];
  
  Double_t gn_bf = 0.5*TMath::Erfc( -Q_bf/( sqrt(2.0)*s_bf ) );
  Double_t k_bf = s_bf/gn_bf/sqrt( 2.0*TMath::Pi() )*TMath::Exp( -pow( Q_bf, 2.0 )/( 2.0*pow( s_bf, 2.0 ) ) );
  Double_t Qg_bf = Q_bf + k_bf;
  Double_t G_bf = ( w_bf/alpha_bf + (1.0-w_bf)*Qg_bf );

  Double_t dG = ( G_bf/G - 1.0 )*100.0;
  
  cout << " True Gain : " << G << endl;
  cout << " BF Gain   : " << G_bf  << endl;
  cout << " Deviation : " << dG << endl;
  
  cout << "" << endl;
  cout << "" << endl;

  if( fit.chi2r<5.0 )
    {
      h_g->Fill( G_bf );
      h_dg->Fill( dG );
      h_chi2->Fill(fit.chi2r);

      i++;
      
    }
    
  if ( i==1 )
    {
      c1->Update();
      c1->WaitPrimitive();
    }

  for ( Int_t i=0; i<npeaks; i++ ) delete grPE[i];

    }

  h_g->SetMarkerStyle( 20 );
  h_g->SetMarkerSize( 0.65 );
  h_g->SetLineColor( kBlack );
  h_g->SetLineWidth( 2.0 );
  h_g->SetMarkerColor( kBlack );
  h_g->Draw( "HIST" );
  
  c1->Update();
  c1->WaitPrimitive();

  h_dg->SetMarkerStyle( 20 );
  h_dg->SetMarkerSize( 0.65 );
  h_dg->SetLineColor( kBlack );
  h_dg->SetLineWidth( 2.0 );
  h_dg->SetMarkerColor( kBlack );
  h_dg->Draw( "HIST" );
  
  c1->Update();
  c1->WaitPrimitive();
  
  h_chi2->SetMarkerStyle( 20 );
  h_chi2->SetMarkerSize( 0.65 );
  h_chi2->SetLineColor( kBlack );
  h_chi2->SetLineWidth( 2.0 );
  h_chi2->SetMarkerColor( kBlack );
  h_chi2->Draw( "HIST" );
  
  c1->Update();
  c1->WaitPrimitive();

  f->cd();
  f->Write();
  f->Close();
  
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
  argc *= 1.0;
  
  TApplication theApp( "App", 0, 0 );
  
  Float_t mu = atof( argv[1] );
  
  Int_t status = mc( mu );

  return status;
  
}
