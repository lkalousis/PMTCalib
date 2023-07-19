
#include "TApplication.h"

#include "iostream"
#include "iomanip"
#include "fstream"

#include "TCanvas.h"

#include "PMTStyle.h"
#include "PMType.h"

#include "Pedestal.h"
#include "SPEResponse.h"
#include "PMT.h"

#include "PMTModel.h"
#include "SPEFitter.h"

#include "TLatex.h"

using namespace std;

void stripCommas(string &a)     
{				
  size_t pos;
  while( (pos = a.find_first_of(",")) != string::npos) a.replace(pos, 1, " ");
  
}

Double_t get_xmax( TString file )
{
  ifstream input;
  string sline; 
  istringstream iss;
  Double_t var1, var2;
  
  input.open( file, ios::in );
  for ( int i=0; i<5; i++ ) getline( input, sline );
  
  getline( input, sline );
  stripCommas( sline );   
  iss.str( sline );                                            
  iss >> var1 >> var2;
  iss.clear();
  
  Double_t xmax = 1.0;
  xmax *= -1.0e+9*var1;

  input.close();
  input.clear();
      
  return xmax; 

}

Double_t get_xmin( TString file )
{
  ifstream input;
  string sline; 
  istringstream iss;
  Double_t var1, var2;
  
  input.open( file, ios::in );
  for ( int i=0; i<5; i++ ) getline( input, sline );
  
  while ( getline( input, sline ) )
  {
    stripCommas( sline );   
    iss.str( sline );                                            
    iss >> var1 >> var2;
    iss.clear();

  }
  
  Double_t xmin = 1.0;
  xmin *= -1.0e+9*var1;

  input.close();
  input.clear();
      
  return xmin; 

}

Double_t get_nbins( TString file )
{
  ifstream input;
  string sline; 
  Int_t nbins = 0;
  
  input.open( file, ios::in );
  for ( int i=0; i<5; i++ ) getline( input, sline );
  
  while ( getline( input, sline ) )
  {
    nbins++;
    
  }
  
  input.close();
  input.clear();
      
  return nbins; 

}

TH1D* get_histo( Int_t nbins, Double_t xmin, Double_t xmax, TString file )
{
  Double_t wbin = ( xmax-xmin )/(1.0*nbins-1.0);
  TH1D *histo = new TH1D( "histo",  "; Charge in nVs; Entries",  nbins, xmin-wbin/2.0, xmax+wbin/2.0 );

  ifstream input;
  string sline; 
  istringstream iss;
  Double_t var1, var2;
  
  input.open( file, ios::in );
  for ( int i=0; i<5; i++ ) getline( input, sline );
  
  while ( getline( input, sline ) )
  {
    stripCommas( sline );   
    iss.str( sline );                                            
    iss >> var1 >> var2;
    iss.clear();

    var1 *= -1.0e+9;
    
    histo->Fill( var1, var2 );

  }
  
  input.close();
  input.clear();

  for ( int i=1; i<=nbins; i++ ) histo->SetBinError( i, sqrt( histo->GetBinContent( i ) ) );
  
  histo->SetMarkerStyle( 20 );
  histo->SetMarkerSize( 0.75 );
  histo->SetLineColor( kBlack );
  histo->SetMarkerColor( kBlack );
  histo->Draw( "PEZ" );
      
  return histo;
  
}

Int_t project8()
{
  time_t start;  
  
  time( &start );
  
  cout << "" << endl;
  cout << " The macro starts ( project11.C ) ... " << endl;
  cout << "" << endl;

  gROOT->Reset();
  PMTStyle::SetDefaultStyle();
    
  TCanvas *c1 = new TCanvas( "c1", "" );
  c1->cd();
  c1->SetLogy();

  
  TH1D *hSG[6];

  TH1F *h_g = new TH1F( "h_g" , "; Gain in nVs; Entries",  40, 0.025, 0.028 );
  TH1F *h_chi2 = new TH1F( "h_chi2" , "; #chi^{2}/NDOF; Entries",  25, 0.0, 5.0 );
  
  for ( Int_t run=0; run<6; run++ )
    {
      cout << " Run : " << run << endl;
      cout << "" << endl;
            
      Int_t nbins;
      Double_t xmax;
      Double_t xmin;
      Double_t wbin;
            
      xmin = get_xmin( Form( "/Users/kalousis/PMTCalib/work/arxiv-2304.07292/data/%d/F1Title00000.txt", run ) );
      xmax = get_xmax( Form( "/Users/kalousis/PMTCalib/work/arxiv-2304.07292/data/%d/F1Title00000.txt", run ) );
      nbins = get_nbins( Form( "/Users/kalousis/PMTCalib/work/arxiv-2304.07292/data/%d/F1Title00000.txt", run ) );

      nbins /= 8; // <---------- !!! Only for even nbins
      wbin = ( xmax-xmin )/(1.0*nbins-1.0);
      cout << xmin << ", " << xmax << ", " << nbins << ", " << wbin << endl;
      //getchar();
      
      hSG[run] = get_histo( nbins, xmin, xmax, Form( "/Users/kalousis/PMTCalib/work/arxiv-2304.07292/data/%d/F1Title00000.txt", run ) );
      hSG[run]->SetName( Form( "hSG_%d", run ) );
      hSG[run]->GetXaxis()->SetRangeUser( -0.04, 0.40 );
      hSG[run]->SetStats( 0 );
      hSG[run]->Draw( "PEZ" );

      //c1->Update();
      //c1->WaitPrimitive();
  
      /* ... */
      
      SPEFitter fit;
      
      Double_t QQ = hSG[run]->GetBinCenter( hSG[run]->GetMaximumBin() );
      Pedestal ped( QQ, 0.002 ); // 0.002 by eye !!! <--------------------------
      ped.LocatePedestal( hSG[run], QQ, 0.002 );
      
      Double_t mu_test = fit.FindMu( hSG[run], ped.Q0, ped.s0 );
      Double_t g_test = fit.FindG( hSG[run], ped.Q0, mu_test );
      Double_t p_test[4] = {  g_test, 0.3*g_test, 1.0/(0.5*g_test), 0.2 }; 
            
      SPEResponse gauss_test( PMType::GAUSS, p_test ); 
      DFTmethod dft( 2.0*nbins, xmin, xmax, gauss_test );
      dft.wbin = hSG[run]->GetBinWidth(1);
      dft.Norm = hSG[run]->Integral();
      dft.Q0 = ped.Q0;
      dft.s0 = ped.s0;
      dft.mu = mu_test;
        
      fit.SetDFTmethod( dft );
      fit.FitwDFTmethod( hSG[run] );
      
      dft.Norm = fit.vals[0];
  
      dft.Q0 = fit.vals[1];
      dft.s0 = fit.vals[2];
      
      dft.mu = fit.vals[3]; 
      
      Double_t p_fit[4] = { fit.vals[4], fit.vals[5], fit.vals[6], fit.vals[7] };
      dft.spef.SetParams( p_fit );

      TGraph *grBF = dft.GetGraph();
      grBF->Draw( "SAME,L" );
      
      /* ... */
      
      Double_t Q = fit.vals[4];
      Double_t s = fit.vals[5];
      Double_t a = fit.vals[6];
      Double_t w = fit.vals[7];
      
      Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );
      Double_t k = s/gn/sqrt( 2.0*TMath::Pi() )*TMath::Exp( -pow( Q, 2.0 )/( 2.0*pow( s, 2.0 ) ) );
      Double_t Qf = Q + k;

      cout << " Q  : " << Q << endl;
      cout << " Qf : " << Qf << ", " <<  ( Q/Qf - 1.0 )*100.0 << endl;
      cout << "" << endl;
      
      Double_t Gfit = w/a + (1.0-w)*Qf;
      cout << " G  : " << Gfit << endl;
      cout << "    : " << ( Q/Gfit - 1.0 )*100.0 << endl;
      cout << "    : " << ( Qf/Gfit - 1.0 )*100.0 << endl;
      cout << "" << endl;
      
      cout << " D. : " << ( Gfit/0.02656 - 1.0 )*100.0 << endl;
      cout << "" << endl;

      TLatex *t = new TLatex();
      t->DrawLatex( 0.20, 2.0e+4, Form( "#font[22]{#scale[0.8]{#mu = %.3f #pm %.3f}}", fit.vals[3], fit.errs[3] ) );
      TLatex *t1 = new TLatex();
      t1->DrawLatex( 0.20, 1.0e+4, Form( "#font[22]{#scale[0.8]{Q = %.5f #pm %.5f}}", fit.vals[4], fit.errs[4] ) );
      TLatex *t2 = new TLatex();
      t2->DrawLatex( 0.20, 5.0e+3, Form( "#font[22]{#scale[0.8]{#sigma = %.5f #pm %.5f}}", fit.vals[5], fit.errs[5] ) );
      TLatex *t3 = new TLatex();
      t3->DrawLatex( 0.20, 2.5e+3, Form( "#font[22]{#scale[0.8]{#alpha = %.0f #pm %.0f}}", fit.vals[6], fit.errs[6] ) );
      TLatex *t4 = new TLatex();
      t4->DrawLatex( 0.20, 1.4e+3, Form( "#font[22]{#scale[0.8]{w = %.3f #pm %.3f}}", fit.vals[7], fit.errs[7] ) );
      
      h_g->Fill( Gfit );
      h_chi2->Fill(fit.chi2r);

      
        Int_t npeaks = 10;
        TGraph *grPE[npeaks];
        
        for ( Int_t i=0; i<npeaks; i++ )
        {
            grPE[i] = dft.GetGraphN( i );
            grPE[i]->Draw( "SAME,L" );
            
        }
        
        //while ( run==4 ) { c1->Update(); c1->WaitPrimitive(); }
        //c1->Update();
        //c1->WaitPrimitive();
        
        for ( Int_t i=0; i<npeaks; i++ ) delete grPE[i];
    }

  h_g->SetMarkerStyle( 20 );
  h_g->SetMarkerSize( 0.65 );
  h_g->SetLineColor( kBlack );
  h_g->SetLineWidth( 2.0 );
  h_g->SetMarkerColor( kBlack );
  h_g->Draw( "HIST" );

  //c1->Update();
  //c1->WaitPrimitive();
  
  h_chi2->SetMarkerStyle( 20 );
  h_chi2->SetMarkerSize( 0.65 );
  h_chi2->SetLineColor( kBlack );
  h_chi2->SetLineWidth( 2.0 );
  h_chi2->SetMarkerColor( kBlack );
  h_chi2->Draw( "HIST" );
  
  //c1->Update();
  //c1->WaitPrimitive();
  
  cout << "" << endl;
  cout << "" << endl;
  
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

int main() 
{
  TApplication theApp( "App", 0, 0 );
  
  Int_t status = project8();

  return status;
  
}


/*

input.close();
      input.clear();

      input.open( Form( "/Users/kalousis/PMTCalib/work/arxiv/data/%d/F1Title00000.txt", run ), ios::in );
      for ( int i=0; i<5; i++ ) getline( input, sline );
      
      while ( getline( input, sline ) ) 
	{
	  stripCommas( sline );   
	  iss.str( sline );                                            
	  
	  iss >> var1 >> var2;
	  iss.clear();
	    
	  var1 *= -1.0e+9;
	  
	  hSG[run]->Fill( var1, var2 );
	  if ( var1<-0.03 ) break;
	  
	  iss.clear();
	  
	}

      input.close();
      input.clear();
           
      for ( int j=1; j<=hSG[run]->GetXaxis()->GetNbins(); j++ ) hSG[run]->SetBinError( j, sqrt( hSG[run]->GetBinContent( j ) ) );
      
      hSG[run]->SetMarkerStyle( 20 );
      hSG[run]->SetMarkerSize( 0.75 );
      hSG[run]->SetLineColor( kBlack );
      hSG[run]->SetMarkerColor( kBlack );
      hSG[run]->Draw( "PEZ" );

      SPEFitter fit;
      
      Double_t QQ = hSG[run]->GetBinCenter( hSG[run]->GetMaximumBin() );
      Pedestal ped( QQ, 0.002 );
      ped.LocatePedestal( hSG[run], QQ, 0.002 );
      
      Double_t mu_test = fit.FindMu( hSG[run], ped.Q0, ped.s0 );
      Double_t g_test = fit.FindG( hSG[run], ped.Q0, mu_test );
      Double_t p_test[4] = {  g_test, 0.3*g_test, 1.0/(0.5*g_test), 0.2 }; // Estimated

      SPEResponse gauss_test( PMType::GAUSS, p_test );
      DFTmethod dft( 2.0*nbins, xmin, xmax, gauss_test );

      dft.wbin = hSG[run]->GetBinWidth(1);
  
      dft.Norm = hSG[run]->Integral();
      dft.Q0 = ped.Q0;
      dft.s0 = ped.s0;
      dft.mu = mu_test;
        
      fit.SetDFTmethod( dft );
      fit.FitwDFTmethod( hSG[run] );
      
      dft.Norm = fit.vals[0];
  
      dft.Q0 = fit.vals[1];
      dft.s0 = fit.vals[2];
      
      dft.mu = fit.vals[3]; 
      
      Double_t p_fit[4] = { fit.vals[4], fit.vals[5], fit.vals[6], fit.vals[7] };
      dft.spef.SetParams( p_fit );

      TGraph *grBF = dft.GetGraph();
      grBF->Draw( "SAME,L" );

      Double_t Q = fit.vals[4];
      Double_t s = fit.vals[5];
      Double_t a = fit.vals[6];
      Double_t w = fit.vals[7];
      
      Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );
      Double_t k = s/gn/sqrt( 2.0*TMath::Pi() )*TMath::Exp( -pow( Q, 2.0 )/( 2.0*pow( s, 2.0 ) ) );
      Double_t Qf = Q + k;

      cout << " Q  : " << Q << endl;
      cout << " Qf : " << Qf << endl;
      cout <<  ( Q/Qf - 1.0 )*100.0 << endl;
      cout << "" << endl;
            
      Double_t Gfit = ( w/a+(1.0-w)*Qf ); 
      cout << " BF Gain : " << Gfit  << endl;
      cout << "" << endl;
      
      h_g->Fill( Gfit );
      h_chi2->Fill(fit.chi2r);
      
      c1->Update();
      c1->WaitPrimitive();
*/
