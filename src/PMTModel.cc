
#include "PMTModel.h"
#include <gsl_sf_hyperg.h>

using namespace std;

ClassImp( PMTModel )

PMTModel::PMTModel()
{}

PMTModel::~PMTModel()
{}

PMTModel::PMTModel( Int_t _nbins, Double_t _xmin, Double_t _xmax, PMType::Model _modtype )
{
  nbins = _nbins;

  xmin = _xmin;
  xmax = _xmax;

  step = ( xmax-xmin )/( 1.0*nbins*1.0 );
  
  modtype = _modtype;
  
  if ( _modtype==PMType::SIMPLEGAUSS ) nparams = 8;
  if ( _modtype==PMType::TRUNCGAUSS ) nparams = 8;
  if ( _modtype==PMType::ANATRUNCG ) nparams = 8;
  
}

void PMTModel::SetParams( Double_t _params[] )
{
  for ( Int_t i=0; i<nparams; i++ )
    {
      params[i] = _params[i];
      
    }
  
  return;
  
}

Double_t PMTModel::GetValue( Double_t xx )
{
  Double_t result = -666;
  
  if ( modtype==PMType::SIMPLEGAUSS ) result = F1( xx );
  if ( modtype==PMType::TRUNCGAUSS ) result = F2( xx );
  if ( modtype==PMType::ANATRUNCG ) result = F3( xx );
  
  return result;
  
}

Double_t PMTModel::F1( Double_t xx )
{
  Double_t result = 0.0; 
  
  Double_t Norm = params[0];
  
  Double_t Q0 = params[1];
  Double_t s0 = params[2];
  
  Double_t mu = params[3];
    
  Double_t Q = params[4];
  Double_t s = params[5];
  
  Double_t alpha = params[6];
  Double_t w = params[7];
    
  Double_t arg0 = 0.0; 
  if ( s0!=0.0 ) arg0 = ( xx - Q0 )/s0;    
  else cout << "Error: The code tries to divide by zero ! " << endl;
  result += TMath::Poisson( 0.0, mu )/( sqrt( 2.0*TMath::Pi() )*s0 )*TMath::Exp( -0.5*arg0*arg0 );

  Double_t Q1 = Q0 + Q;
  Double_t s1 = sqrt( pow( s0, 2.0 ) + pow( s, 2.0 ) );
  Double_t arg1 = 0.0; 
  if ( s1!=0.0 ) arg1 = ( xx - Q1 )/s1;    
  else cout << "Error: The code tries to divide by zero ! " << endl;
  
  Double_t omega = ( Q0+alpha*pow( s0, 2.0 )-xx )/sqrt(2.0)/s0;
  Double_t SR1 = w*alpha/2.0*TMath::Exp( -alpha*( xx-Q0 )+pow( alpha*s0, 2.0 )/2.0 )*TMath::Erfc( omega );
  SR1 += (1.0-w)/( sqrt( 2.0*TMath::Pi() )*s1 )*TMath::Exp( -0.5*arg1*arg1 );
  result += TMath::Poisson( 1.0, mu )*SR1;

  for ( Int_t n=2; n<25; n++ )
    {
      Double_t Qs = w/alpha+(1.0-w)*Q;
      Double_t Qn = Q0 + 1.0*n*Qs;
      Double_t ss2 = w/pow( alpha, 2.0 ) + (1-w)*pow( s, 2.0 ) + w*(1.0-w)*pow( Q-1.0/alpha, 2.0 );
      Double_t sn2 = pow( s0, 2.0 ) + 1.0*n*ss2;
      Double_t sn = sqrt( sn2 );
                  
      Double_t argn = 0.0; 
      if ( sn!=0.0 ) argn = ( xx-Qn )/sn;    
      else cout << "Error: The code tries to divide by zero ! " << endl;
      result += TMath::Poisson( 1.0*n, mu )/( sqrt( 2.0*TMath::Pi() )*sn ) * TMath::Exp( -0.5*argn*argn );
                  
    }
  
  result *= Norm*wbin;
  
  return result;

}

Double_t PMTModel::F2( Double_t xx )
{
  Double_t result = 0.0; 
  
  Double_t Norm = params[0];
  
  Double_t Q0 = params[1];
  Double_t s0 = params[2];
  
  Double_t mu = params[3];
    
  Double_t Q = params[4];
  Double_t s = params[5];
  
  Double_t alpha = params[6];
  Double_t w = params[7];
    
  Double_t arg0 = 0.0; 
  if ( s0!=0.0 ) arg0 = ( xx - Q0 )/s0;    
  else cout << "Error: The code tries to divide by zero ! " << endl;
  result += TMath::Poisson( 0.0, mu )/( sqrt( 2.0*TMath::Pi() )*s0 )*TMath::Exp( -0.5*arg0*arg0 );

  Double_t Q1 = Q0 + Q;
  Double_t s1 = sqrt( pow( s0, 2.0 ) + pow( s, 2.0 ) );
  Double_t arg1 = 0.0; 
  if ( s1!=0.0 ) arg1 = ( xx - Q1 )/s1;    
  else cout << "Error: The code tries to divide by zero ! " << endl;

  Double_t omega = ( Q0+alpha*pow( s0, 2.0 )-xx )/sqrt(2.0)/s0;
  Double_t SR1 = w*alpha/2.0*TMath::Exp( -alpha*( xx-Q0 )+pow( alpha*s0, 2.0 )/2.0 )*TMath::Erfc( omega );
  
  Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );
  Double_t A = ( Q0-xx )*pow( s, 2.0 ) - Q*pow( s0, 2.0 ); 
  Double_t B = sqrt( 2.0 )*s0*s*s1;
  SR1 += (1.0-w)*1.0/2.0/gn/( sqrt( 2.0*TMath::Pi() )*s1 )*TMath::Exp( -0.5*arg1*arg1 )*TMath::Erfc( A/B );
  result += TMath::Poisson( 1.0, mu )*SR1;

  Double_t k = s/gn/sqrt( 2.0*TMath::Pi() )*TMath::Exp( -pow( Q, 2.0 )/( 2.0*pow( s, 2.0 ) ) );
  Double_t Qf = Q + k;
  Double_t sf2 = pow( s, 2.0 ) - ( Q+k )*k;
        
  Double_t Qs = w/alpha+(1.0-w)*Qf;
  Double_t ss2 = w/pow( alpha, 2.0 ) + (1-w)*sf2 + w*(1.0-w)*pow( Qf-1.0/alpha, 2.0 );
  
  for ( Int_t n=2; n<25; n++ )
    {
      Double_t Qn = Q0 + 1.0*n*Qs;
      Double_t sn2 = pow( s0, 2.0 ) + 1.0*n*ss2;
      Double_t sn = sqrt( sn2 );
                  
      Double_t argn = 0.0; 
      if ( sn!=0.0 ) argn = ( xx-Qn )/sn;    
      else cout << "Error: The code tries to divide by zero ! " << endl;
      result += TMath::Poisson( 1.0*n, mu )/( sqrt( 2.0*TMath::Pi() )*sn ) * TMath::Exp( -0.5*argn*argn );
            
    }
  
  result *= Norm*wbin;
  
  return result;
  
}

Double_t PMTModel::F3( Double_t xx )
{
  Double_t result = 0.0; 
  
  Double_t Norm = params[0];
  
  Double_t Q0 = params[1];
  Double_t s0 = params[2];
  
  Double_t mu = params[3];
    
  Double_t Q = params[4];
  Double_t s = params[5];
  
  Double_t alpha = params[6];
  Double_t w = params[7];

  /* ... */

  Double_t arg = 0.0; 
  if ( s0!=0.0 ) arg = ( xx - Q0 )/s0;    
  else cout << "Error: The code tries to divide by zero." << endl;
    
  Double_t SR0 = 1.0/( sqrt( 2.0*TMath::Pi() )*s0 )*TMath::Exp( -0.5*arg*arg );
  SR0 *= TMath::Poisson( 0, mu );
  result += SR0; // 0

  
  Double_t omega0 = ( xx - Q0 - alpha*pow( s0, 2.0 ) )/sqrt(2.0)/s0;
  Double_t SR1 = w*alpha/2.0*TMath::Exp( -alpha*( xx-Q0 )+pow( alpha*s0, 2.0 )/2.0 )*TMath::Erfc( -omega0 );

  Double_t Q1 = Q0+Q;
  Double_t s12 = pow( s0, 2.0 )+pow( s, 2.0 );
  Double_t s1 = sqrt( s12 );
  
  Double_t arg1 = 0.0; 
  if ( s1!=0.0 ) arg1 = ( xx - Q1 )/s1;    
  else cout << "Error: The code tries to divide by zero." << endl;

  Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );
  Double_t A = ( Q0-xx )*pow( s, 2.0 ) - Q*pow( s0, 2.0 ); 
  Double_t B = sqrt( 2.0 )*s0*s*s1;
  SR1 += ( 1.0-w )/2.0/gn/( sqrt( 2.0*TMath::Pi() )*s1 )*TMath::Exp( -0.5*arg1*arg1 )*TMath::Erfc( A/B );
  SR1 *= TMath::Poisson( 1, mu );
  result += SR1; // 1
  
  
  Int_t nlim = 10;

  Double_t k = s/gn/sqrt( 2.0*TMath::Pi() )*TMath::Exp( -pow( Q, 2.0 )/( 2.0*pow( s, 2.0 ) ) );
  Double_t Qg = Q + k;
  Double_t sg2 = pow( s, 2.0 ) - ( Q+k )*k;
  
  for ( Int_t n = 2; n<nlim; n++ )
    {
      Double_t SRn = 0.0;
            
      Double_t Qn = Q0 + 1.0*n*Qg;
      Double_t sn2 = pow( s0,2.0 )+ 1.0*n*sg2;
      Double_t sn = sqrt( sn2 );
      
      Double_t argn = 0.0; 
      if ( sn!=0.0 ) argn = ( xx - Qn )/sn;    
      else cout << "Error: The code tries to divide by zero." << endl;
      Double_t gnB = 1.0/( sqrt( 2.0*TMath::Pi() )*sn )*TMath::Exp( -0.5*argn*argn );
      SRn += pow( 1.0-w, n )*gnB;
      
      
      for ( Int_t m=1; m<=n; m++ )
	{
	  Double_t Qmn = Q0 + 1.0*(n-m)*Qg;
	  Double_t smn2 = pow( s0, 2.0 )+1.0*(n-m)*sg2;
	  Double_t smn = sqrt( smn2 );
	  
	  Double_t cmn = alpha*pow( alpha*smn*sqrt( 2.0 ), m-1.0 )/TMath::Factorial( m-1.0 )/2.0/sqrt( TMath::Pi() );
	  
	  Double_t psi = ( xx-Qmn )/sqrt(2.0)/smn;
	  Double_t psi2 = pow( psi, 2.0 );
	  Double_t omega = ( xx-Qmn-alpha*pow( smn, 2.0 ) )/sqrt(2.0)/smn;
	  Double_t omega2 = pow( omega, 2.0 );
	  
	  Double_t A1m = 1.0*m/2.0;
	  Double_t A2m = (0.0+m+1.0)/2.0;
	  
	  Double_t Imn=0.0;
	  Double_t hi_limit=25.0;
	  
	  if ( omega>=hi_limit )
	    {
	      Imn = 2.0*sqrt( TMath::Pi() )*TMath::Exp( omega2-psi2 + ( m-1.0 )*TMath::Log( omega ) );
	      
	    }
	  else if ( omega<hi_limit && omega>=0.0  )
	    {
	      Double_t t1 = TMath::Gamma( A1m )*gsl_sf_hyperg_1F1( 1.0/2.0-A1m, 1.0/2.0, -omega2 );
	      Double_t t2 = 2.0*omega*TMath::Gamma( A2m )*gsl_sf_hyperg_1F1( 3.0/2.0-A2m, 3.0/2.0, -omega2 );
	      Imn = ( t1+t2 )*TMath::Exp( omega2-psi2 );
	      
	    }
	  else if ( omega<0.0 )
	    {
	      Double_t t3 = TMath::Gamma( A1m )*TMath::Gamma( A2m )/sqrt( TMath::Pi() );
	      Imn = t3*gsl_sf_hyperg_U( A1m, 1.0/2.0, omega2 )*TMath::Exp( -psi2 );
	      
	    }
	  
	  
	  Double_t hmnB = cmn*Imn;
	  Double_t binom = TMath::Factorial( n )/TMath::Factorial( m )/TMath::Factorial( n-m );
	  SRn += binom*pow( w, m )*pow( 1.0-w, n-m )*hmnB;
	  
	}
           
      SRn *= TMath::Poisson( n, mu );
      result += SRn; // n= 2-nlim
      
    } 
  
  
  Double_t Qs = w/alpha + (1.0-w)*Qg;
  Double_t ss2 = w/pow( alpha, 2.0 ) + (1-w)*sg2 + w*(1.0-w)*pow( Qg-1.0/alpha, 2.0 );
      
  for ( Int_t n = nlim; n<50; n++ )
    {
      Double_t Qn = Q0 + 1.0*n*Qs;
      Double_t sn2 = pow( s0, 2.0 ) + 1.0*n*ss2;
      Double_t sn = sqrt( sn2 );
      
      Double_t argn = 0.0; 
      if ( sn!=0.0 ) argn = ( xx - Qn )/sn;    
      else cout << "Error: The code tries to divide by zero." << endl;
      Double_t SRn = 1.0/( sqrt( 2.0*TMath::Pi() )*sn )*TMath::Exp( -0.5*argn*argn );

      SRn *= TMath::Poisson( n, mu );
      result += SRn; // n >= nlim
      
    }
      
  

  /* ... */
    
  result *= Norm*wbin;
  
  return result;
  
}

TGraph* PMTModel::GetGraph()
{
  Double_t x[nbins];
  Double_t y[nbins];
  
  for ( Int_t i=0; i<nbins; i++ )
    {
      x[i] = xmin + step/2 + 1.0*i*step;

      Double_t y_ = GetValue( x[i] );
      
      if ( y_<1.0e-10 ) y[i] = 1.e-3;
      else y[i] = y_;
      
    }
  
  TGraph *_gr = new TGraph( nbins, x, y );
  
  _gr->SetLineWidth( 2 );
  int cc = kAzure+1;
  _gr->SetLineColor( cc );
  _gr->SetMarkerColor( cc );
  _gr->SetMarkerSize( 0.1 );
  
  return _gr;

}


  
      




/*

//Double_t cn = alpha*pow( alpha*s0*sqrt( 2.0 ), n-1.0 )/TMath::Factorial( n-1.0 )/2.0/sqrt( TMath::Pi() );

      Double_t psi = ( xx-Q0 )/sqrt(2.0)/s0;
      Double_t psi2 = pow( psi, 2.0 );
      Double_t omega02 = pow( omega0, 2.0 );

      Double_t A1 = 1.0*n/2.0;
      Double_t A2 = (0.0+n+1.0)/2.0;
      
      Double_t I0=1.0;
      Double_t hi_limit = 20.0;
      //Double_t lo_limit = -1000.0;
      
      if ( omega02>=hi_limit )
	{
	  I0 = 2.0*sqrt( TMath::Pi() )*TMath::Exp( omega02-psi2 + ( n-1.0 )*TMath::Log( omega0 ) );
	  
	}
      else if ( omega02<hi_limit && omega02>=0.0  )
	{
	  Double_t t1 = TMath::Gamma( A1 )*gsl_sf_hyperg_1F1( 1.0/2.0-A1, 1.0/2.0, -omega02 );
	  Double_t t2 = 2.0*omega0*TMath::Gamma( A2 )*gsl_sf_hyperg_1F1( 3.0/2.0-A2, 3.0/2.0, -omega02 );
	  I0 = ( t1+t2 )*TMath::Exp( omega02-psi2 );
	  
	}
      else if ( omega02<0.0 )
	{
	  Double_t t3 = TMath::Gamma( A1 )*TMath::Gamma( A2 )/sqrt( TMath::Pi() );
	  I0 = t3*gsl_sf_hyperg_U( A1, 1.0/2.0, omega02 )*TMath::Exp( -psi2 );
      
	}

*/
