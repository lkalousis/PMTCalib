
#include "TApplication.h"

#include "iostream"
#include "iomanip"

#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLine.h"
#include "TGraph.h"

#include "PMTStyle.h"
#include "PMType.h"

#include "Pedestal.h"
#include "SPEResponse.h"
#include "PMT.h"

#include "DFTmethod.h"

using namespace std;

Int_t plot2( Float_t mu )
{
  time_t start;  
  time( &start );
  cout << "" << endl;
  cout << " The macro starts ( plot2.C ) ... " << endl;
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

  Int_t nbins = 2000;
  Double_t xmin = -0.14;
  Double_t xmax = 4000.0*0.0005+xmin;
  
  PMT specimen( nbins, xmin, xmax, ped, gamma );
  Int_t ntot = 2.5e+5;
  
  specimen.GenSpectrum( ntot, mu );
  specimen.GetSpectrum()->SetStats(0);
  //specimen.DrawSpectrum();
  
  //c1->Update();
  //c1->WaitPrimitive();
  
  DFTmethod dft( 2.0*nbins, xmin, xmax, gamma );
      
  dft.wbin = specimen.GetSpectrum()->GetBinWidth(1);
  dft.Norm = ntot;
  
  dft.Q0 = ped.Q0;
  dft.s0 = ped.s0;
      
  dft.mu = mu;

  TGraph *grBF = dft.GetGraph();
  grBF->SetLineColor( kBlack );
  grBF->GetXaxis()->SetRangeUser( -0.05, 0.11 );
  grBF->SetMinimum(0.0);
  grBF->SetMaximum(25000);
  grBF->SetTitle( "Pedestal; Charge in A.U.; " );
  grBF->Draw( "AL" );
  
  //c1->Update();
  //c1->WaitPrimitive();
  

  TH1D *h = new TH1D( "h", "", nbins, xmin, xmax );
  Double_t sum = 0.0;

  for ( Int_t i=1; i<=nbins; i++ )
    {
      Double_t x1 = h->GetXaxis()->GetBinCenter( i );
      //Double_t x2 = h->GetXaxis()->GetBinUpEdge( i );
      Double_t y1 = grBF->Eval( x1 );
      sum += y1; 
      
    }

  cout << sum << endl;
  
  Double_t x[3000];
  Double_t y[3000];
  Int_t n = 0;

  x[0] = h->GetXaxis()->GetBinLowEdge( 1 );
  y[0] = 0.0;
  n++;

  Double_t sum1 = 0.0;
  
  for ( Int_t i=1; i<=nbins; i++ )
    {
      Double_t x1 = h->GetXaxis()->GetBinCenter( i );
      Double_t x2 = h->GetXaxis()->GetBinUpEdge( i );
      Double_t y1 = grBF->Eval( x1 );
      sum1 += y1/sum; 

      x[i] = x2;
      y[i] = sum1;
      n++;
      
    }

  TGraph *gr1 = new TGraph( n, x, y );
  gr1->SetLineWidth( 2 );
  //gr1->SetLineColor( kRed );
  gr1->GetXaxis()->SetRangeUser( -0.05, 0.11 );
  gr1->SetMinimum(0.0);
  gr1->SetMaximum(1.0);
  gr1->SetTitle( "Cumulative distribution; Charge in A.U.; " );
  gr1->Draw("AL");

  TGraph *gr2 = new TGraph( n, y, x );
  Double_t cut = gr2->Eval( 0.10 );
  cout << cut << endl;
    
  TLine *l1 = new TLine( cut, 0.0, cut, 0.1 );
  l1->SetLineWidth( 2 );
  l1->SetLineColor( kRed );
  l1->SetLineStyle(9);
  l1->Draw();

  TLine *l2 = new TLine( -0.05, 0.1, cut, 0.1 );
  l2->SetLineWidth( 2 );
  l2->SetLineColor( kRed );
  l2->SetLineStyle(9);
  l2->Draw();

  TLine *l3 = new TLine( Q0+5.0*s0, 0.0, Q0+5.0*s0, 1.0 );
  l3->SetLineWidth( 2 );
  l3->SetLineColor( kBlack );
  l3->SetLineStyle(9);
  l3->Draw();
  
  Double_t x5[1];
  Double_t y5[1];
  x5[0] = cut;
  y5[0] = 0.1;

  TGraph *gr5 = new TGraph( 1, x5, y5 );
  gr5->SetMarkerColor( kRed );
  gr5->SetMarkerStyle(20);
  gr5->SetMarkerSize(1);
  gr5->Draw( "P,SAME" );
  
  c1->Update();
  c1->WaitPrimitive();
  

  grBF->Draw( "AL" );

  TLine *l4 = new TLine( Q0+5.0*s0, 0.0, Q0+5.0*s0, 25000.0 );
  l4->SetLineWidth( 2 );
  //l4->SetLineColor(kRed);
  l4->SetLineStyle(9);
  l4->Draw();

  TLine *l5 = new TLine( cut, 0.0, cut, 25000.0 );
  l5->SetLineWidth( 2 );
  l5->SetLineColor(kRed);
  l5->SetLineStyle(9);
  l5->Draw();

  
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

int main( int argc, char *argv[] )  
{
  argc++;
  Float_t mu = atof( argv[1] );
    
  TApplication theApp( "App", 0, 0 );
  Int_t status = plot2( mu );
  return status;
  
}

