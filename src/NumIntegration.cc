
#include "NumIntegration.h"

using namespace std;


ClassImp( NumIntegration )

NumIntegration::NumIntegration()
{}

NumIntegration::~NumIntegration()
{}

NumIntegration::NumIntegration( Int_t _nbins, Double_t _xmin, Double_t _xmax, SPEResponse _spef )
{
  nbins = _nbins;
  
  xmin = _xmin;
  xmax = _xmax;

  step = ( xmax-xmin )/( 1.0*nbins*1.0 );
  
  spef = _spef;
  
  N = nbins; 
  
  xvalues.clear();
  
  for ( UInt_t i=0; i<N; i++ )
    {
      Double_t xx = xmin + wbin/2.0 + 1.0*i*step;
      
      xvalues.push_back( xx );
      
    }
  
}

void NumIntegration::CalculateValues()
{
  Double_t x[N];
  Double_t y[N];
  
  yvalues.clear();
  
  for ( UInt_t i=0; i<N; i++ )
    {
      Double_t xx = xvalues.at( i );
      x[i] = xx;
      
      Double_t result = 0.0;

      Double_t arg = 0.0; 
      if ( s0!=0.0 ) arg = ( xx - Q0 )/s0;    
      else cout << "Error: The code tries to divide by zero (num1) " << endl;
      Double_t SR0 = TMath::Exp(-mu)/( sqrt( 2.0*TMath::Pi() )*s0 ) * TMath::Exp( -0.5*arg*arg );
      result += SR0;

      Double_t x_lo = Q0-7.0*s0;
      Double_t x_hi = Q0+7.0*s0;
      UInt_t nds = 100;
      Double_t ds = ( x_hi-x_lo )/nds;
      
      Double_t SR1 = 0.0;
	
      for ( UInt_t j=0; j<nds; j++ )
	{
	  Double_t tt = x_lo + j*ds;
	  Double_t xt = xx - tt;

	  Double_t arg = 0.0; 
	  if ( s0!=0.0 ) arg = ( tt - Q0 )/s0;    
	  else cout << "Error: The code tries to divide by zero (num2) " << endl;
	  Double_t yy0 = 1.0/( sqrt( 2.0*TMath::Pi() )*s0 ) * TMath::Exp( -0.5*arg*arg );
	  	  
	  if ( xt>=0.0 ) SR1 += yy0*spef.GetValue( xt )*ds;
	  
	}

      SR1 *= TMath::Poisson( 1, mu );
      result += SR1;

      Double_t Q = spef.params[0];
      Double_t s = spef.params[1];
      Double_t alpha = spef.params[2];
      Double_t w = spef.params[3];
      
      Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );
      Double_t k = s/gn/sqrt( 2.0*TMath::Pi() )*TMath::Exp( -pow( Q, 2.0 )/( 2.0*pow( s, 2.0 ) ) );
      Double_t Qg = Q + k;
      Double_t sg2 = pow( s, 2.0 ) - ( Q+k )*k;
      
      Double_t Qs = w/alpha + (1.0-w)*Qg;
      Double_t ss2 = w/pow( alpha, 2.0 ) + (1-w)*sg2 + w*(1.0-w)*pow( Qg-1.0/alpha, 2.0 );
      
      for ( Int_t n = 2; n<50; n++ )
	{
	  Double_t Qn = Q0 + 1.0*n*Qs;
	  Double_t sn2 = pow( s0, 2.0 ) + 1.0*n*ss2;
	  Double_t sn = sqrt( sn2 );
	  
	  Double_t argn = 0.0; 
	  if ( sn!=0.0 ) argn = ( xx - Qn )/sn;    
	  else cout << "Error: The code tries to divide by zero." << endl;
	  Double_t SRn = 1.0/( sqrt( 2.0*TMath::Pi() )*sn )*TMath::Exp( -0.5*argn*argn );
	  
	  SRn *= TMath::Poisson( n, mu );
	  result += SRn; 
	  
	}
       
      result *= wbin*Norm;
      
      if ( isnan(result) || isinf(result) )
	{
	  cout << "yep" << endl;
	  //getchar();
	  
	}
	  
      yvalues.push_back( result );
      y[i]=result;
      
    }

  gr = new TGraph( N, x, y );
  
  return;
  
}

Double_t NumIntegration::GetValue( Double_t xx )
{
  Double_t y_ = gr->Eval( xx ); 

  return y_;

}

TGraph* NumIntegration::GetGraph()
{
  CalculateValues();
  
  Double_t x[N];
  Double_t y[N];
  
  for ( UInt_t i=0; i<N; i++ )
    {
      x[i] = xvalues.at( i );

      Double_t y_ = GetValue( x[i] );
      
      if ( y_<1.0e-10 ) y[i] = 1.e-4;
      else y[i] = y_;
      
    }
  
  TGraph *_gr = new TGraph( N, x, y );
  
  _gr->SetLineWidth( 2 );
  int cc = kAzure+1;
  _gr->SetLineColor( cc );
  _gr->SetMarkerColor( cc );
  _gr->SetMarkerSize( 0.1 );
    
  return _gr;
  
}


/*
  Int_t nlim = 2; // !!!
  TGraph *gr_s[nlim];
  
  Double_t x[N];
  Double_t y_s[nlim][N];
    
  for ( UInt_t i=0; i<N; i++ )
    {
      x[i] = xvalues.at( i );
      
      Double_t arg = 0.0; 
      if ( s0!=0.0 ) arg = ( x[i] - Q0 )/s0;    
      else cout << "Error: The code tries to divide by zero (num) " << endl;
      Double_t yy = 1.0/( sqrt( 2.0*TMath::Pi() )*s0 ) * TMath::Exp( -0.5*arg*arg );
      y_s[0][i] = yy;
      
      y_s[1][i] = spef.GetValue( x[i] );
      
    }
  
  gr_s[0] = new TGraph( N, x, y_s[0] );
  gr_s[1] = new TGraph( N, x, y_s[1] );
  
  
  
  for ( Int_t k=2; k<nlim; k++ )
    {
      //Int_t k=2;
      //cout << k << endl;
      for ( UInt_t i=0; i<N; i++ )
	{
	  Double_t xx = xvalues.at( i );
	  Double_t result = 0.0;
          
	  for ( UInt_t j=0; j<=i; j++ )
	    {
	      Double_t tt = xvalues.at( j );
	      Double_t xt = xx - tt;
	      if ( xt>=0.0 ) result += spef.GetValue( tt )*gr_s[k-1]->Eval( xt )*step;
	      //result += gr_s[1]->Eval( tt )*gr_s[k-1]->Eval( xx-tt )*step; // this must be wrong !!!
	      	  
	    }
	  
	  y_s[k][i] = result;
	  
	}
      
      gr_s[k] = new TGraph( N, x, y_s[k] );
      
    }
    
  Double_t y_SID[N];
  
  for ( UInt_t i=0; i<N; i++ )
    {
      Double_t xx = xvalues.at( i );
      y_SID[i] = 0.0;
      
      for ( Int_t k=1; k<nlim; k++ )
	{
	  if ( xx>=0.0 ) y_SID[i] += TMath::Poisson( k, mu )*gr_s[k]->Eval( xx );
	  
	}
  
    }

  
  
  
  TGraph *gr_SID = new TGraph( N, x, y_SID );
  
  
  Double_t y_SR[N];
  
  for ( UInt_t i=0; i<N; i++ )
    {
      Double_t xx = xvalues.at( i );
      Double_t result = 0.0;

       Double_t x_lo = Q0-5.0*s0;
      Double_t x_hi = Q0+5.5*s0;
      UInt_t M = ( x_hi-x_lo )/step;
      
      //for ( UInt_t j=0; j<=i; j++ )
      for ( UInt_t j=0; j<=M; j++ )
	{
	  //Double_t tt = xvalues.at( j );
	  Double_t tt = x_lo + j*step;
	  Double_t xt = xx - tt;
	  
	  Double_t arg = 0.0; 
	  if ( s0!=0.0 ) arg = ( tt - Q0 )/s0;    
	  else cout << "Error: The code tries to divide by zero (num) " << endl;
	  Double_t yy0 = 1.0/( sqrt( 2.0*TMath::Pi() )*s0 ) * TMath::Exp( -0.5*arg*arg );
	  	  
	  if ( xt>=0.0 ) result += yy0*gr_SID->Eval( xt )*step;
	  //result += yy0*gr_SID->Eval( xx-tt )*step;
	  
	}
      
      y_SR[i] = result;
      
    }
  
  TGraph *gr_SR = new TGraph( N, x, y_SR );
  
  yvalues.clear();
  
  Double_t Q = spef.params[0];
  Double_t s = spef.params[1];
  Double_t alpha = spef.params[2];
  Double_t w = spef.params[3];
  
  Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );
  Double_t k = s/gn/sqrt( 2.0*TMath::Pi() )*TMath::Exp( -pow( Q, 2.0 )/( 2.0*pow( s, 2.0 ) ) );
  Double_t Qf = Q + k;
  Double_t sf2 = pow( s, 2.0 ) - ( Q+k )*k;

  Double_t Qs = w/alpha+(1.0-w)*Qf;
  Double_t ss2 = w/pow( alpha, 2.0 ) + (1-w)*sf2 + w*(1.0-w)*pow( Qf-1.0/alpha, 2.0 );
  
  for ( UInt_t i=0; i<N; i++ )
    {
      Double_t xx = xvalues.at( i );

      Double_t arg = 0.0; 
      if ( s0!=0.0 ) arg = ( xx - Q0 )/s0;    
      else cout << "Error: The code tries to divide by zero (num) " << endl;
      Double_t yy0 = 1.0/( sqrt( 2.0*TMath::Pi() )*s0 ) * TMath::Exp( -0.5*arg*arg );

      Double_t yy = TMath::Exp(-mu)*( yy0 ) +  gr_SR->Eval( xx );
      
      for ( UInt_t j=nlim; j<75; j++ )
	{
	  Double_t Qn = Q0 + 1.0*j*Qs;
	  Double_t sn2 = pow( s0, 2.0 ) + 1.0*j*ss2;
	  Double_t sn = sqrt( sn2 );
                  
	  Double_t argn = 0.0; 
	  if ( sn!=0.0 ) argn = ( xx-Qn )/sn;    
	  else cout << "Error: The code tries to divide by zero ! " << endl;
	  yy += TMath::Poisson( 1.0*j, mu )/( sqrt( 2.0*TMath::Pi() )*sn ) * TMath::Exp( -0.5*argn*argn );
	  	  
	}
      
      yy *= Norm*wbin;
      
      yvalues.push_back( yy );
      
    }
  
  Double_t y[N];

  for ( UInt_t i=0; i<N; i++ )
    {
      y[i] = yvalues.at( i );
      
    }

  gr = new TGraph( N, x, y );
  
  for ( Int_t k=0; k<nlim; k++ ) delete gr_s[k];
  delete gr_SID;
  delete gr_SR;


 */
