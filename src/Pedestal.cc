
#include "Pedestal.h"

Double_t wpin0;

using namespace std;

ClassImp( Pedestal )

Double_t _pedfunc( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];
    
  Double_t Q_0 = par[0];
  Double_t s_0 = par[1];
        
  Double_t arg = 0.0; 
  if ( s_0!=0.0 ) arg = ( xx - Q_0 )/s_0;    
  else cout << "Error: The code tries to divide by zero." << endl;
  
  Double_t result = 1.0/( sqrt( 2.0*TMath::Pi() ) * s_0 ) * TMath::Exp( -0.5*arg*arg );
    
  return result;
  
}

Double_t _pedfit( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];

  Double_t N = par[0];
  
  Double_t Q_0 = par[1];
  Double_t s_0 = par[2];
        
  Double_t arg = 0.0; 
  if ( s_0!=0.0 ) arg = ( xx - Q_0 )/s_0;    
  //else cout << "Error: The code tries to divide by zero." << endl;
  
  Double_t result = wpin0*N/( sqrt( 2.0*TMath::Pi() )*s_0 ) * TMath::Exp( -0.5*arg*arg );
    
  return result;
  
}

Pedestal::Pedestal()
{}

Pedestal::~Pedestal()
{}

Pedestal::Pedestal( Double_t _Q0, Double_t _s0 )
{
  Q0 = _Q0;
  s0 = _s0;
  
  pedfunc = new TF1( "pedfunc", _pedfunc, Q0-25.0*s0, Q0+25.0*s0, 2 );
  pedfunc->SetLineColor( kRed );
  pedfunc->SetLineWidth( 2.0 );
  pedfunc->SetNpx( 10000 );
  pedfunc->SetParameters( Q0, s0 );
  
}

Double_t Pedestal::GenQ()
{
  Double_t _x = pedfunc->GetRandom();

  return _x;
  
}

void Pedestal::LocatePedestal( TH1 *hspec, Double_t _Q0, Double_t _s0 )
{
  pedfit = new TF1( "pedfit", _pedfit, _Q0-10.0*_s0, _Q0-10.0*_s0, 3 );
  pedfit->SetParNames( "Norm", "Q", "#sigma" );
  pedfit->SetLineColor( kRed );
  pedfit->SetNpx( 10000 );
  
  wpin0 = hspec->GetXaxis()->GetBinWidth(1);
  
  Double_t Norm = hspec->Integral();
  
  pedfit->SetParameter( 0, Norm );
  pedfit->SetParLimits( 0, Norm*0.005, Norm*500.0 );

  pedfit->SetParameter( 1, _Q0 );
  pedfit->SetParLimits( 1, _Q0-5.0*_s0, _Q0+5.0*_s0 );

  pedfit->SetParameter( 2, _s0 );
  pedfit->SetParLimits( 2, _s0*0.1, _s0*10.0 );
  
  status = hspec->Fit( "pedfit", "", "", _Q0-3.0*_s0, _Q0+2.6*_s0 );
  //status = hspec->Fit( "pedfit", "Q0", "", _Q0-8.0*_s0, _Q0+8.0*_s0 ); // Config. or MUTEL 
  chi2 = pedfit->GetChisquare();
  //pedfit->GetNDF();
  
  Double_t *pars0 = pedfit->GetParameters();
  const Double_t *erpars0 = pedfit->GetParErrors();
  
  Q0 = pars0[1];
  s0 = pars0[2];
  
  dQ0 = erpars0[1];
  ds0 = erpars0[2];
  
  return;
  
}
