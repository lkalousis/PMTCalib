
#include "TApplication.h"

#include "iostream"
#include "iomanip"
#include "fstream"
#include "string"
#include "sstream"

#include "TString.h"
#include "TLatex.h"
#include "PMTStyle.h"
#include "SPEFitter.h"
#include "Pedestal.h"
#include "SPEResponse.h"
#include "DFTmethod.h"
#include "Occupancy.h"

#include "TH1D.h"
#include "TCanvas.h"


using namespace std;

TString file;
Int_t nbins;
Double_t x1;
Double_t xmin;
Double_t x2;
Double_t xmax;
Double_t wbin0;
Double_t wbin;
istringstream iss;
ifstream input;
string sline; 

Double_t var1, var2;

void stripCommas(string &a)     
{				
  size_t pos;
  while( (pos = a.find_first_of(",")) != string::npos) a.replace(pos, 1, " ");
  
}

void get_x1()
{
  input.open( file, ios::in );
  for ( int i=0; i<5; i++ ) getline( input, sline );
  
  while ( getline( input, sline ) )
  {
    stripCommas( sline );   
    iss.str( sline );                                            
    iss >> var1 >> var2;
    iss.clear();

  }
  
  x1 = 1.0;
  x1 *= -1.0e+9*var1;

  input.close();
  input.clear();
      
  return;

}

void get_x2()
{
  input.open( file, ios::in );
  for ( int i=0; i<5; i++ ) getline( input, sline );
  
  getline( input, sline );
  stripCommas( sline );   
  iss.str( sline );                                            
  iss >> var1 >> var2;
  iss.clear();
  
  x2 = 1.0;
  x2 *= -1.0e+9*var1;

  input.close();
  input.clear();
  
  return;

}

void get_nbins()
{  
  input.open( file, ios::in );
  for ( int i=0; i<5; i++ ) getline( input, sline );

  nbins = 0;
  
  while ( getline( input, sline ) )
  {
    nbins++;
    
  }
  
  input.close();
  input.clear();
      
  return; 

}

TH1D* get_histo()
{
  nbins /= 2;
  TH1D *histo = new TH1D( "_h",  "; Charge in nVs; Entries",  nbins, xmin, xmax );
    
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
  //histo->Draw( "PEZ" );
      
  return histo;
  
}

Int_t ana()
{
  time_t start;  
  
  time( &start );
  
  cout << "" << endl;
  cout << " The macro starts ( ana.C ) ... " << endl;
  cout << "" << endl;
  
  //gROOT->Reset();
  PMTStyle::SetDefaultStyle();
  
  TCanvas *c1 = new TCanvas( "c1", "" );
  c1->cd();
  c1->SetLogy();

  Int_t run[15] = { 10, 14, 18, 20, 24, 28, 30, 34, 38, 40, 44, 48, 50, 54, 58 };
  
  TH1D *hSG[15];
  
  TH1F *h_g = new TH1F( "h_g" , "; Gain in nVs; Entries",  50, 0.085, 0.105 );
  TH1F *h_g2 = new TH1F( "h_g2" , "; Gain in nVs; Entries",  50, 0.085, 0.105 );
  
  TH1F *h_chi2 = new TH1F( "h_chi2" , "; #chi^{2}/NDOF; Entries",  25, 0.0, 5.0 );
    
  for ( Int_t i=0; i<15; i++ )
    {
      cout << " i   : " << i << endl;
      cout << " run : " << run[i] << endl;
      
      file = Form( "/Users/kalousis/PMTCalib/work/arxiv-2607.24250/data-iv/%d/single/F1Title00000.txt", run[i] );
      cout << file << endl;

      get_x1();
      cout << x1 << endl;
      
      get_x2();
      cout << x2 << endl;

      get_nbins();
      cout << nbins << endl;

      wbin0 = (x2-x1)/(nbins-1.0);
      cout << wbin0 << endl;

      xmin = x1-wbin0/2.0;
      xmax = x2+wbin0/2.0;

      hSG[i] = get_histo();
      hSG[i]->SetName( Form( "hSG_%d", i ) );
      hSG[i]->GetXaxis()->SetRangeUser( -0.05, 1.25 );
      hSG[i]->SetStats( 0 );
      hSG[i]->Draw( "PEZ" );

      wbin = hSG[i]->GetXaxis()->GetBinWidth(1);
      cout << wbin << endl;
      cout << hSG[i]->GetXaxis()->GetBinLowEdge(1) << endl;
      cout << hSG[i]->GetXaxis()->GetBinUpEdge(nbins/2) << endl;

      Int_t iQ0 = hSG[i]->GetMaximumBin();
      Double_t max = hSG[i]->GetBinContent( iQ0 );
      hSG[i]->SetMinimum(0.1);
      hSG[i]->SetMaximum( max*2.5 );
      
      //c1->Update();
      //c1->WaitPrimitive();

      /* ... */

      SPEFitter fit;

      Double_t Q0 = hSG[i]->GetBinCenter( iQ0 );
      Pedestal ped( Q0, 0.002 );
      ped.LocatePedestal( hSG[i], Q0, 0.002 );
      cout << ped.Q0 << ", " << ped.s0 << endl;
      
      Double_t mu_test = fit.FindMu( hSG[i], ped.Q0, ped.s0 );
      Double_t g_test = fit.FindG( hSG[i], ped.Q0, mu_test );
      cout << mu_test << ", " << g_test << endl;

      Double_t p_test[4] = { 1.0/g_test, 7.0, 1.0/(0.5*g_test), 0.3 };
      SPEResponse gamma_test( PMType::GAMMA, p_test );

      DFTmethod dft( nbins, xmin, xmax, gamma_test );
      
      dft.wbin = hSG[i]->GetBinWidth(1);
      dft.Norm = hSG[i]->Integral();
  
      dft.Q0 = ped.Q0;
      dft.s0 = ped.s0;
      
      dft.mu = mu_test;
        
      fit.SetDFTmethod( dft );
      fit.FitwDFTmethod( hSG[i] );

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
      
  
      Float_t Gfit = ( fit.vals[7]/fit.vals[6]+(1.0-fit.vals[7])/fit.vals[4] ); 
      cout << " BF Gain : " << Gfit  << endl;
      //cout << "" << endl;

      if ( fit.chi2r<3.0 )
	{
	  h_g->Fill( Gfit );
	  h_chi2->Fill(fit.chi2r);

	  Occupancy occ( fit.vals[1], fit.vals[2] );
	  Float_t Gfit2 = occ.FindG( hSG[i], 0.2 );
	  h_g2->Fill( Gfit2 );
	   
	  cout << " BF Gain : " << Gfit2  << endl;
	  cout << "" << endl;
	  
	}
      
      cout << Form( "%.3f", fit.vals[3] ) << " $\\pm$ " <<  Form( "%.3f", fit.errs[3] ) << " & "
	   << Form( "%.2f", fit.vals[7] ) << " $\\pm$ " <<  Form( "%.2f", fit.errs[7] ) << " & "
	   << Form( "%.0f", fit.vals[6] ) << " $\\pm$ " <<  Form( "%.0f", fit.errs[6] ) << " & "
	   << Form( "%.2f", fit.vals[4] ) << " $\\pm$ " <<  Form( "%.2f", fit.errs[4] ) << " & "
	   << Form( "%.1f", fit.vals[5] ) << " $\\pm$ " <<  Form( "%.1f", fit.errs[5] ) << " & "
	   << Form( "%.2f", fit.chi2r   ) << " \\\\ "  
	   << endl;
      
      cout << "" << endl;

      c1->Update();
      c1->WaitPrimitive();

      delete grBF;
      for ( Int_t i=0; i<npeaks; i++ ) delete grPE[i];
      
      /* ... */

      cout << "" << endl;
      cout << "" << endl;
      
    }

  h_g2->SetMarkerStyle( 20 );
  h_g2->SetMarkerSize( 0.65 );
  h_g2->SetLineColor( kBlack );
  h_g2->SetLineWidth( 2.0 );
  h_g2->SetMarkerColor( kBlack );
  h_g2->Draw( "HIST" );

  c1->Update();
  c1->WaitPrimitive();
  
  h_g->SetMarkerStyle( 20 );
  h_g->SetMarkerSize( 0.65 );
  h_g->SetLineColor( 50 );
  h_g->SetLineWidth( 2.0 );
  h_g->SetMarkerColor( 50 );
  h_g->Draw( "SAMES,HIST" );
  
  while (1!=0){
  c1->Update();
  c1->WaitPrimitive();
  }
  
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
