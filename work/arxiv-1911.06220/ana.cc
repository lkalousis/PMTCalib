
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

Int_t ana()
{
  time_t start;  
  
  time( &start );
  
  cout << "" << endl;
  cout << " The macro starts ( ana.C ) ... " << endl;
  cout << "" << endl;

  gROOT->Reset();
  PMTStyle::SetDefaultStyle();
    
  TCanvas *c1 = new TCanvas( "c1", "" );
  c1->cd();
  c1->SetLogy();

  
  TH1D *hSG[61];

  TH1F *h_g = new TH1F( "h_g" , "; Gain in nVs; Entries",  40, 0.092, 0.097 );
  TH1F *h_chi2 = new TH1F( "h_chi2" , "; #chi^{2}/NDOF; Entries",  25, 0.0, 5.0 );
  
  for ( Int_t run=0; run<61; run++ )
    {
      cout << " Run : " << run << endl;
      cout << "" << endl;
            
      Int_t nbins;
      Double_t xmax;
      Double_t xmin;
      Double_t wbin;

      ifstream input;
      input.open( Form( "/Users/kalousis/PMTCalib/work/arxiv-1911.06220/data/%d/single/F1Title00000.txt", run ), ios::in );
      if ( ( !input ) || ( input.peek()==EOF )  )	      
	{
	  input.close();
	  input.clear();
	  cout << " Warning : input file does not exist or it is empty" << endl;
	  cout << "" << endl;
	  continue;
	  
	}
      
      xmin = get_xmin( Form( "/Users/kalousis/PMTCalib/work/arxiv-1911.06220/data/%d/single/F1Title00000.txt", run ) );
      xmax = get_xmax( Form( "/Users/kalousis/PMTCalib/work/arxiv-1911.06220/data/%d/single/F1Title00000.txt", run ) );
      nbins = get_nbins( Form( "/Users/kalousis/PMTCalib/work/arxiv-1911.06220/data/%d/single/F1Title00000.txt", run ) );

      nbins /= 2; // <---------- !!! Only for even nbins
      wbin = ( xmax-xmin )/(1.0*nbins-1.0);
      //cout << xmin << ", " << xmax << ", " << nbins << ", " << wbin << endl;
      //getchar();
      
      hSG[run] = get_histo( nbins, xmin, xmax, Form( "/Users/kalousis/PMTCalib/work/arxiv-1911.06220/data/%d/single/F1Title00000.txt", run ) );
      hSG[run]->SetName( Form( "hSG_%d", run ) );
      hSG[run]->GetXaxis()->SetRangeUser( -0.1, 2.0 );
      hSG[run]->SetStats( 0 );
      hSG[run]->Draw( "PEZ" );

      /* ... */

      
      SPEFitter fit;
      
      Double_t qqq = hSG[run]->GetBinCenter( hSG[run]->GetMaximumBin() );
      Pedestal ped( qqq, 0.002 );
      ped.LocatePedestal( hSG[run], qqq, 0.002 );
      
      
      Double_t mu_test = fit.FindMu( hSG[run], ped.Q0, ped.s0 );
      Double_t g_test = fit.FindG( hSG[run], ped.Q0, mu_test );
      
      Double_t p_test[4] = { 1.0/g_test, 5.0, 1.0/(0.1*g_test), 0.2 };
      SPEResponse gaus_test( PMType::GAMMA, p_test );
      
      DFTmethod dft( 2.0*nbins, xmin, xmax, gaus_test );
      
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

      Int_t npeaks = 10;
      TGraph *grPE[npeaks];
      
      for ( Int_t i=0; i<npeaks; i++ )
	{
	  grPE[i] = dft.GetGraphN( i );
	  grPE[i]->Draw( "SAME,L" );
	  
	}
      
      Double_t Gfit = ( fit.vals[7]/fit.vals[6]+(1.0-fit.vals[7])/fit.vals[4] ); 
      cout << " BF Gain : " << Gfit  << endl;
      cout << "" << endl;

      if ( fit.chi2r<3.0 )
	{
	  h_g->Fill( Gfit );
	  h_chi2->Fill(fit.chi2r);

	}
      
      /* ... */
            
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
  
  Int_t status = ana();

  return status;
  
}


/*
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
	
*/
      
//while ( run==4 ) { c1->Update(); c1->WaitPrimitive(); } 
