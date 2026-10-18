
#include "SPEResponse.h"

using namespace std;

ClassImp( SPEResponse )


Double_t _gausexpfunc( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];
    
  Double_t Q = par[0];
  Double_t s = par[1];

  Double_t alpha = par[2];
  Double_t w = par[3];
    
  Double_t result = 0.0;
  if ( xx>=0.0 )
    {
      Double_t arg = 0.0; 
      if ( s!=0.0 ) arg = ( xx - Q )/s;    
      else cout << "Error: The code tries to divide by zero." << endl;

      Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );
      
      result = w*alpha*TMath::Exp( -xx*alpha ) + ( 1.0-w )/( sqrt( 2.0*TMath::Pi() )*s*gn )*TMath::Exp( -0.5*arg*arg );

    }
  
  return result;
  
}

Double_t _gaus2( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];
    
  Double_t Q1 = par[0];
  Double_t s1 = par[1];
    
  //Double_t alpha = par[2];
  //Double_t w1 = par[3];
  
  //Double_t w2 = par[4];
  //Double_t Q2 = par[5];
  //Double_t s2 = par[6];

  Double_t Q2 = par[2];
  Double_t s2 = par[4];
  Double_t w = par[3];
  
  Double_t result = 0.0;
  if ( xx>=0.0 )
    {
      Double_t arg1 = 0.0; 
      if ( s1!=0.0 ) arg1 = ( xx - Q1 )/s1;    
      else cout << "Error: The code tries to divide by zero." << endl;

      Double_t arg2 = 0.0; 
      if ( s2!=0.0 ) arg2 = ( xx - Q2 )/s2;    
      else cout << "Error: The code tries to divide by zero." << endl;
      
      Double_t gn1 = 0.5*TMath::Erfc( -Q1/( sqrt(2.0)*s1 ) );
      Double_t gn2 = 0.5*TMath::Erfc( -Q2/( sqrt(2.0)*s2 ) );
      
      result = //w1*alpha*TMath::Exp( -xx*alpha ) +
	w/( sqrt( 2.0*TMath::Pi() )*s2*gn2 )*TMath::Exp( -0.5*arg2*arg2 )
	+ ( 1.0-w )/( sqrt( 2.0*TMath::Pi() )*s1*gn1 )*TMath::Exp( -0.5*arg1*arg1 );

    }
  
  return result;
  
}

Double_t _gaus2expfunc( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];
    
  Double_t Q = par[0];
  Double_t s = par[1];

  Double_t alpha1 = par[2];
  Double_t w1 = par[3];
  Double_t alpha2 = par[4];
  Double_t w2 = par[5];
  
  Double_t result = 0.0;
  if ( xx>=0.0 )
    {
      Double_t arg = 0.0; 
      if ( s!=0.0 ) arg = ( xx - Q )/s;    
      else cout << "Error: The code tries to divide by zero." << endl;

      Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );
      
      result = w1*alpha1*TMath::Exp( -xx*alpha1 ) + w2*alpha2*TMath::Exp( -xx*alpha2 ) + ( 1.0-w1-w2 )/( sqrt( 2.0*TMath::Pi() )*s*gn )*TMath::Exp( -0.5*arg*arg );

    }
  
  return result;
  
}

Double_t _gammaexpfunc( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];
    
  Double_t lambda = par[0];
  Double_t theta = par[1];
  
  Double_t alpha = par[2];
  Double_t w = par[3];
    
  Double_t result = 0.0;
  if ( xx>=0.0 )
    {
      Double_t f = lambda*(1.0+theta);
      Double_t fx = lambda*(1.0+theta)*xx;
      
      result = w*1.0*alpha*TMath::Exp( -xx*alpha ) + ( 1.0-w )*f*pow( fx, theta )/TMath::Gamma(1.0+theta)*TMath::Exp( -fx );

    }
  
  return result;
  
}

Double_t _gamma2expfunc( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];
    
  Double_t lambda = par[0];
  Double_t theta = par[1];
  
  Double_t alpha1 = par[2];
  Double_t w1 = par[3];

  Double_t alpha2 = par[4];
  Double_t w2 = par[5];
    
  Double_t result = 0.0;
  if ( xx>=0.0 )
    {
      Double_t f = lambda*(1.0+theta);
      Double_t fx = lambda*(1.0+theta)*xx;
      
      result = w1*1.0*alpha1*TMath::Exp( -xx*alpha1 )+w2*1.0*alpha2*TMath::Exp( -xx*alpha2 )+( 1.0-w1-w2 )*f*pow( fx, theta )/TMath::Gamma(1.0+theta)*TMath::Exp( -fx );

    }
  
  return result;
  
}

Double_t _weibullexpfunc( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];
  
  Double_t lambda = par[0];
  Double_t kappa = par[1];
  
  Double_t alpha = par[2];
  Double_t w = par[3];
    
  Double_t result = 0.0;
  if ( xx>=0.0 )
    {
      Double_t f = kappa/lambda;
      Double_t fx = xx/lambda;
      
      result = w*alpha*TMath::Exp( -xx*alpha ) + ( 1.0-w )*f*pow( fx, kappa-1.0 )*TMath::Exp( -pow( fx, kappa ) );
      
    }
  
  return result;
  
}

Double_t _lognormalexpfunc( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];
    
  Double_t Q = par[0];
  Double_t s = par[1];

  Double_t alpha = par[2];
  Double_t w = par[3];
    
  Double_t result = 0.0;
  if ( xx>0.0 )
    {
      Double_t arg = 0.0; 
      if ( s!=0.0 ) arg = ( TMath::Log( xx ) - Q )/s;    
      else cout << "Error: The code tries to divide by zero." << endl;
            
      result = w*alpha*TMath::Exp( -xx*alpha ) + ( 1.0-w )/( sqrt( 2.0*TMath::Pi() )*s*xx )*TMath::Exp( -0.5*arg*arg );

    }
  
  return result;
  
}

Double_t _testfunc( Double_t *x, Double_t *par )
{
   Double_t xx = x[0];
    
  Double_t Q = par[0];
  Double_t s = par[1];
  
  Double_t lambda = par[2];
  Double_t theta = par[4];
  Double_t w = par[3];
    
  Double_t result = 0.0;
  if ( xx>=0.0 )
    {
      Double_t arg = 0.0; 
      if ( s!=0.0 ) arg = ( xx - Q )/s;    
      else cout << "Error: The code tries to divide by zero." << endl;

      Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );

      Double_t f = lambda*(1.0+theta);
      Double_t fx = lambda*(1.0+theta)*xx;
            
      result = w*f*pow( fx, theta )/TMath::Gamma(1.0+theta)*TMath::Exp( -fx ) + ( 1.0-w )/( sqrt( 2.0*TMath::Pi() )*s*gn )*TMath::Exp( -0.5*arg*arg );

    }
  
  return result;
  
}

SPEResponse::SPEResponse()
{}

SPEResponse::~SPEResponse()
{}

SPEResponse::SPEResponse( PMType::Response _spetype, Double_t _params[] )
{
  spetype = _spetype;
    
  
  if ( spetype==PMType::GAUSS )
    {
      params[0] = _params[0];
      params[1] = _params[1];
      
      params[2] = _params[2];
      params[3] = _params[3];
      
      spefunc = new TF1( "spefunc", _gausexpfunc, params[0]-80.0*params[1], params[0]+80.0*params[1], 4 );
      spefunc->SetParameters( params[0], params[1], params[2], params[3] );

      nparams = 4;
      
    }
  
  if ( spetype==PMType::GAMMA )
    {
      params[0] = _params[0];
      params[1] = _params[1];
      
      params[2] = _params[2];
      params[3] = _params[3];
      
      spefunc = new TF1( "spefunc", _gammaexpfunc, 1.0/params[0]-80.0*1.0/params[0]/sqrt(1.0+params[1]), 1.0/params[0]+80.0*1.0/params[0]/sqrt(1.0+params[1]), 4 );
      spefunc->SetParameters( params[0], params[1], params[2], params[3] );

      nparams = 4;
      
    }

  if ( spetype==PMType::WEIBULL )
    {
      params[0] = _params[0];
      params[1] = _params[1];
      
      params[2] = _params[2];
      params[3] = _params[3];
      
      spefunc = new TF1( "spefunc", _weibullexpfunc, params[0]-80.0*params[0], params[0]+80.0*params[0], 4 );
      spefunc->SetParameters( params[0], params[1], params[2], params[3] );

      nparams = 4;
      
    }


  if ( spetype==PMType::LOGNORMAL )
    {
      params[0] = _params[0];
      params[1] = _params[1];
      
      params[2] = _params[2];
      params[3] = _params[3];
      
      spefunc = new TF1( "spefunc", _lognormalexpfunc, TMath::Log( params[0] )-80.0*TMath::Log( params[0] ), TMath::Log( params[0] )+80.0*TMath::Log( params[0] ), 4 );
      spefunc->SetParameters( params[0], params[1], params[2], params[3] );

      nparams = 4;
      
    }

   if ( spetype==PMType::GAUSS2EXP )
    {
      params[0] = _params[0];
      params[1] = _params[1];
      
      params[2] = _params[2];
      params[3] = _params[3];

      params[4] = _params[4];
      params[5] = _params[5];
      
      spefunc = new TF1( "spefunc", _gaus2expfunc, params[0]-80.0*params[1], params[0]+80.0*params[1], 6 );
      spefunc->SetParameters( params[0], params[1], params[2], params[3], params[4], params[5] );

      nparams = 6;
      
    }

   if ( spetype==PMType::GAMMA2EXP )
    {
      params[0] = _params[0];
      params[1] = _params[1];
      
      params[2] = _params[2];
      params[3] = _params[3];

      params[4] = _params[4];
      params[5] = _params[5];
      
      spefunc = new TF1( "spefunc", _gamma2expfunc, 1.0/params[0]-80.0*1.0/params[0]/sqrt(1.0+params[1]), 1.0/params[0]+80.0*1.0/params[0]/sqrt(1.0+params[1]), 6 );
      spefunc->SetParameters( params[0], params[1], params[2], params[3], params[4], params[5] );

      nparams = 6;
      
    }

    if ( spetype==PMType::GAUSS2 )
    {
      params[0] = _params[0];
      params[1] = _params[1];
      
      params[2] = _params[2];
      params[3] = _params[3];
      params[4] = _params[4];
      //params[5] = _params[5];
      //params[6] = _params[6];
            
      spefunc = new TF1( "spefunc", _gaus2, params[0]-80.0*params[1], params[0]+80.0*params[1], 5 );
      spefunc->SetParameters( params[0], params[1], params[2], params[3], params[4] );//, params[5], params[6] );

      nparams = 5;
      
    }
    
   if ( spetype==PMType::TEST )
    {
      params[0] = _params[0];
      params[1] = _params[1];
      
      params[2] = _params[2];
      params[3] = _params[3];
      params[4] = _params[4];
      
      spefunc = new TF1( "spefunc", _testfunc, params[0]-80.0*params[1], params[0]+80.0*params[1], 5 );
      spefunc->SetParameters( params[0], params[1], params[2], params[3], params[4] );

      nparams = 5;
      
    }
  
    
  spefunc->SetLineColor( kBlue );
  spefunc->SetLineWidth( 2.0 );
  spefunc->SetNpx( 10000 );
  
}

void SPEResponse::SetParams( Double_t _params[] )
{
  for ( Int_t i=0; i<nparams; i++ )
    {
      params[i] = _params[i];
      spefunc->SetParameter( i, params[i] );
      
    }
      
}

Double_t SPEResponse::GetValue( Double_t xx )
{
  Double_t result = spefunc->Eval( xx );

  return result;
  
}


Double_t SPEResponse::GenQ()
{
  Double_t _x = spefunc->GetRandom();

  return _x;
  
}
